// Re-derive the KKT verifier's dual sign violation per column/row and classify it.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include "Highs.h"
static bool hasL(double v){return v>-kHighsInf;} static bool hasU(double v){return v<kHighsInf;}
static double signViol(double lo,double up,double val,double adj,double nearTol){
  bool l=hasL(lo),u=hasU(up);
  if(!l&&!u) return std::fabs(adj);
  if(l&&!u) return std::max(0.0,-adj);
  if(!l&&u) return std::max(0.0,adj);
  if(up-lo<=nearTol) return 0.0;
  double mid=0.5*lo+0.5*up;
  return val<mid?std::max(0.0,-adj):std::max(0.0,adj);
}
int main(int argc,char**argv){
  Highs h; h.setOptionValue("output_flag",false); h.readModel(argv[1]);
  HighsLp lp=h.getLp(); lp.a_matrix_.ensureColwise();
  FILE*f=std::fopen(argv[2],"r"); int nc,nr; if(fscanf(f,"%d %d",&nc,&nr)!=2) return 1;
  std::vector<double> x(nc),d(nc),y(nr);
  for(auto&v:x) if(fscanf(f,"%lf",&v)!=1) return 1;
  for(auto&v:d) if(fscanf(f,"%lf",&v)!=1) return 1;
  for(auto&v:y) if(fscanf(f,"%lf",&v)!=1) return 1;
  double sense = lp.sense_==ObjSense::kMinimize?1.0:-1.0, tol=1e-4, nearTol=std::sqrt(tol);
  std::vector<double> rc(nc), act(nr,0.0);
  for(int j=0;j<nc;++j){ double s=0; for(int e=lp.a_matrix_.start_[j];e<lp.a_matrix_.start_[j+1];++e){ s+=lp.a_matrix_.value_[e]*y[lp.a_matrix_.index_[e]]; act[lp.a_matrix_.index_[e]]+=lp.a_matrix_.value_[e]*x[j]; } rc[j]=lp.col_cost_[j]-s; }
  struct R{double v;int idx;bool col;};
  std::vector<R> all;
  std::map<std::string,std::pair<int,double>> byType; // count with viol>1e-4, sum of viol
  auto type=[&](double lo,double up){ bool l=hasL(lo),u=hasU(up); if(l&&u) return std::string(up-lo<=nearTol?"fixed":"boxed"); if(l) return std::string("lower-only"); if(u) return std::string("upper-only"); return std::string("free");};
  for(int j=0;j<nc;++j){ double v=signViol(lp.col_lower_[j],lp.col_upper_[j],x[j],sense*rc[j],nearTol); all.push_back({v,j,true}); if(v>1e-4){auto&t=byType["col "+type(lp.col_lower_[j],lp.col_upper_[j])]; t.first++; t.second+=v;} }
  for(int i=0;i<nr;++i){ double v=signViol(lp.row_lower_[i],lp.row_upper_[i],act[i],sense*y[i],nearTol); all.push_back({v,i,false}); if(v>1e-4){auto&t=byType["row "+type(lp.row_lower_[i],lp.row_upper_[i])]; t.first++; t.second+=v;} }
  std::sort(all.begin(),all.end(),[](const R&a,const R&b){return a.v>b.v;});
  std::printf("entries with sign violation > 1e-4 by type: \n");
  for(auto&kv:byType) std::printf("   %-16s count %6d  sum viol %.3e\n",kv.first.c_str(),kv.second.first,kv.second.second);
  std::printf("top violations (the verifier divides by 1+near-zero-cost-norm):\n");
  for(int k=0;k<8;++k){ auto&r=all[k]; if(r.col){ int j=r.idx; std::printf("  col %6d viol %.3e  bounds [%g, %g] x=%.4g  rc=%.4g  col_dual(reported)=%.4g cost=%g\n",j,r.v,lp.col_lower_[j],lp.col_upper_[j],x[j],rc[j],d[j],lp.col_cost_[j]); } else { int i=r.idx; std::printf("  row %6d viol %.3e  bounds [%g, %g] activity=%.6g y=%.4g\n",i,r.v,lp.row_lower_[i],lp.row_upper_[i],act[i],y[i]); } }
  return 0;
}
