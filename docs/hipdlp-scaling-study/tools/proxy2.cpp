// Per-(instance, scaling): sigma(A'), ||x'*-x'0||, ||y'*||, omega, K_bal = sig*Dx*Dy, K_w = sig*(w Dx^2 + Dy^2/w)/2
// usage: proxy2 MPS SOL NAME   (x0 = projection of 0 onto bounds, y0 = 0, omega = (1+||c||)/(1+||b||) unscaled)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
using namespace hipdlp::model;
static double l2(const std::vector<double>& v){ long double s=0; for(double x:v) s+=(long double)x*x; return std::sqrt((double)s);}
static SparseMatrix transposeOf(const SparseMatrix& A){ SparseMatrix T; T.num_rows=A.num_cols; T.num_cols=A.num_rows; T.starts.assign(T.num_rows+1,0);
  for(int e=0;e<A.numNonzeros();++e) T.starts[A.indices[e]+1]++; for(int i=0;i<T.num_rows;++i) T.starts[i+1]+=T.starts[i];
  T.indices.resize(A.numNonzeros()); T.values.resize(A.numNonzeros()); std::vector<int> pos(T.starts.begin(),T.starts.end()-1);
  for(int r=0;r<A.num_rows;++r) for(int e=A.starts[r];e<A.starts[r+1];++e){int p=pos[A.indices[e]]++; T.indices[p]=r; T.values[p]=A.values[e];} return T;}
static double sigmaMax(const SparseMatrix& A,const SparseMatrix& AT,int iters=60){ if(A.numNonzeros()==0) return 0; std::vector<double> x(A.num_cols),y(A.num_rows);
  for(int i=0;i<A.num_cols;++i) x[i]=1.0+0.37*std::sin(1.0*i); double n=l2(x),est=0; for(auto&v:x) v/=n;
  for(int it=0;it<iters;++it){ for(int r=0;r<A.num_rows;++r){double s=0; for(int e=A.starts[r];e<A.starts[r+1];++e) s+=A.values[e]*x[A.indices[e]]; y[r]=s;}
    for(int c=0;c<A.num_cols;++c){double s=0; for(int e=AT.starts[c];e<AT.starts[c+1];++e) s+=AT.values[e]*y[AT.indices[e]]; x[c]=s;}
    n=l2(x); est=std::sqrt(n); if(n==0) break; for(auto&v:x) v/=n;} return est;}
static double clamp0(double lo,double up){ return std::min(std::max(0.0,lo),up); }
int main(int argc,char**argv){
  if(argc<4) return 2;
  Highs h; h.setOptionValue("output_flag",false); if(h.readModel(argv[1])==HighsStatus::kError) return 3; const HighsLp lp=h.getLp();
  FILE*f=std::fopen(argv[2],"r"); if(!f) return 4; int nc,nr; if(fscanf(f,"%d %d",&nc,&nr)!=2) return 1;
  std::vector<double> x(nc),d(nc),y(nr); for(auto&v:x) if(fscanf(f,"%lf",&v)!=1) return 1; for(auto&v:d) if(fscanf(f,"%lf",&v)!=1) return 1; for(auto&v:y) if(fscanf(f,"%lf",&v)!=1) return 1;
  ProblemBuildResult built=buildProblem(lp); if(!built.ok()) return 5; const Problem& p=built.problem; const auto& tr=p.transformation();
  HighsLp lc=lp; lc.a_matrix_.ensureColwise(); std::vector<double> act(nr,0.0);
  for(int j=0;j<nc;++j) for(int e=lc.a_matrix_.start_[j];e<lc.a_matrix_.start_[j+1];++e) act[lc.a_matrix_.index_[e]]+=lc.a_matrix_.value_[e]*x[j];
  const int m=p.numRows(), n=p.numCols();
  std::vector<double> xs(n,0.0), x0(n,0.0), ys(m,0.0);
  for(int j=0;j<nc;++j) xs[j]=x[j];
  for(int i=0;i<nr;++i){ int s=tr.slack_col_for_original_row[i]; if(s>=0) xs[s]=act[i]; ys[tr.original_to_transformed_row[i]]=y[i]; }
  for(int j=0;j<n;++j) x0[j]=clamp0(p.colLower()[j],p.colUpper()[j]);
  const double omega=(1.0+p.unscaledNorms().objective_l2)/(1.0+p.unscaledNorms().row_bounds_l2);
  struct Cfg{const char*name;ScalingMode mode;int mask;};
  const Cfg cfgs[]={{"O",ScalingMode::kOff,0},{"C",ScalingMode::kCurtisReid,0},{"A",ScalingMode::kAuto,5},{"AR",ScalingMode::kAuto,1},{"AP",ScalingMode::kAuto,4},
    {"L",ScalingMode::kLegacy,5},{"R",ScalingMode::kLegacy,1},{"P",ScalingMode::kLegacy,4}};
  for(const Cfg&c:cfgs){
    std::vector<double> R(m,1.0),C(n,1.0);
    if(c.mode!=ScalingMode::kOff){ ScalingOptions o; o.mode=c.mode; o.legacy_scaling_mode=c.mask; ScalingResult r=computeScalingPipeline(p.matrix(),o); if(!r.ok()) continue; R=r.scaling.row_scale; C=r.scaling.col_scale; }
    SparseMatrix A=p.matrix(); for(int r=0;r<m;++r) for(int e=A.starts[r];e<A.starts[r+1];++e) A.values[e]*=R[r]*C[A.indices[e]];
    double sig=sigmaMax(A,transposeOf(A));
    std::vector<double> dx(n),dy(m); for(int j=0;j<n;++j) dx[j]=(xs[j]-x0[j])/C[j]; for(int i=0;i<m;++i) dy[i]=ys[i]/R[i];
    double Dx=l2(dx),Dy=l2(dy);
    std::printf("%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",argv[3],c.name,sig,Dx,Dy,omega,sig*Dx*Dy,sig*0.5*(omega*Dx*Dx+Dy*Dy/omega));
  }
  return 0;
}
