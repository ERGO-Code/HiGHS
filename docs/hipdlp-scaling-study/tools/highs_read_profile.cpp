#include <chrono>
#include <cstdio>
#include <ctime>
#include <memory>
#include "Highs.h"
#include "io/Filereader.h"
int main(int argc,char** argv) {
  if(argc!=2) return 2;
  using Clock=std::chrono::steady_clock;
  Highs highs;
  highs.setOptionValue("output_flag",false);
  HighsModel model;
  auto start=Clock::now();
  auto cpu_start=std::clock();
  std::unique_ptr<Filereader> reader(Filereader::getFilereader(highs.getOptions().log_options,argv[1]));
  if(!reader) return 3;
  const auto ret=reader->readModelFromFile(highs.getOptions(),argv[1],model);
  reader.reset();
  if(ret!=FilereaderRetcode::kOk && ret!=FilereaderRetcode::kWarning) return 4;
  double parse=std::chrono::duration<double>(Clock::now()-start).count();
  double parse_cpu=double(std::clock()-cpu_start)/CLOCKS_PER_SEC;
  start=Clock::now();cpu_start=std::clock();
  if(highs.passModel(std::move(model))==HighsStatus::kError) return 5;
  double assess=std::chrono::duration<double>(Clock::now()-start).count();
  double assess_cpu=double(std::clock()-cpu_start)/CLOCKS_PER_SEC;
  std::printf("{\"parse_and_matrix_build_wall_sec\":%.17g,\"parse_cpu_sec\":%.17g,\"pass_model_validation_wall_sec\":%.17g,\"validation_cpu_sec\":%.17g,\"rows\":%lld,\"cols\":%lld,\"nnz\":%zu}\n",parse,parse_cpu,assess,assess_cpu,(long long)highs.getLp().num_row_,(long long)highs.getLp().num_col_,highs.getLp().a_matrix_.value_.size());
}
