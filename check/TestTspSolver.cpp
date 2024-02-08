#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsCallback.h"

#include <fstream>
#include <sstream>

const bool kDebugReport = false;
const HighsInt kmaxIter = 100;

class TspData {
 public:
  TspData() {}
  std::string name_;
  std::string type_;
  HighsInt dimension_;
  std::vector<std::vector<HighsInt>> i_j_to_k_;
  std::vector<std::pair<HighsInt, HighsInt>> k_to_i_j_;
  std::vector<std::vector<HighsInt>> tours_;
  void initialise(std::string& filename, HighsLp& lp);
  void getTours(const double* solution);
  void addCuts(HighsLp& lp);

};

const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

HighsCallbackFunctionType userDefineLazyConstraints =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      TspData tsp_data = *(static_cast<TspData*>(user_callback_data));
      tsp_data.getTours(data_out->mip_solution);
      if (dev_run) {
	printf("TSP: %s has dimension %d and solution has %d tours\n",
	       tsp_data.name_.c_str(),
	       int(tsp_data.dimension_),
	       int(tsp_data.tours_.size()));
      }
	assert(343 == 545);
    };

TEST_CASE("tsp-p01", "[highs_test_tsp_solver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/p01.tsp";
  TspData tsp_data;
  HighsLp lp;
  tsp_data.initialise(filename, lp); 
  
  const double optimal_obective_value = 263;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  //  highs.readModel(filename);
  highs.passModel(lp);
  //  MipData user_callback_data;
  TspData* p_user_callback_data = &tsp_data;
  highs.setCallback(userDefineLazyConstraints, p_user_callback_data);
  printf("Calling highs.setCallback\n");
  highs.startCallback(kCallbackMipDefineLazyConstraints);
  highs.run();
  REQUIRE(highs.getObjectiveValue() == optimal_obective_value);
}

void TspData::initialise(std::string& filename, HighsLp& lp) {
  std::string word;
  std::string format = "";
  HighsInt expected_num_distances;
  std::vector<HighsInt> distances;
  std::ifstream f;
  f.open(filename.c_str(), std::ios::in);
  for (;;) {
    std::stringstream ss;
    std::string strline;
    std::string key;
    bool eof = false;
    if (getline(f, strline)) {
      ss.str(std::string());
      ss << strline; 
      ss >> key;
    } else {
      exit(1);
    }
    if (key == "NAME:" || key == "NAME") {
      if (key == "NAME") {
	ss >> word;
	assert(word == ":");
      }
      ss >> name_;
    } else if (key == "TYPE:" || key == "TYPE") {
      if (key == "TYPE") {
	ss >> word;
	assert(word == ":");
      }
      ss >> type_;
    } else if (key == "DIMENSION:" || key == "DIMENSION") {
      if (key == "DIMENSION") {
	ss >> word;
	assert(word == ":");
      }
      ss >> dimension_;
    } else if (key == "EDGE_WEIGHT_FORMAT:" || key == "EDGE_WEIGHT_FORMAT") {
      if (key == "EDGE_WEIGHT_FORMAT") {
	ss >> word;
	assert(word == ":");
      }
      ss >> format;
    } else if (key == "EDGE_WEIGHT_SECTION") {
      if (format == "") {
	printf("Distance format is not defined\n");
	exit(1);
      } else if (format == "LOWER_DIAG_ROW") {
	expected_num_distances = (dimension_ * (dimension_+1)) / 2;
      } else if (format == "FULL_MATRIX") {
	expected_num_distances = dimension_ * dimension_;
      } else {
	printf("Unknown distance format: %s\n", format.c_str());
	exit(1);
      }
      HighsInt line = 0;
      for (;;) {
	std::stringstream wt_ss;
	std::string wt_strline;
	if (getline(f, wt_strline)) {
	  wt_ss.str(std::string());
	  wt_ss << wt_strline; 
	  HighsInt v;
	  for (;;) {
	    wt_ss >> v;
	    if (wt_ss.tellp() == -1) {
	      line++;
	      break;
	    }
	    distances.push_back(v);
	  }
	  if (HighsInt(distances.size()) == expected_num_distances) break;
	} else {
	  break;
	}
      }
      eof = true;
    }
    if (eof) break;
  }
  printf("TSP: %s\n", name_.c_str());
  printf("TYPE: %s\n", type_.c_str());
  printf("DIMENSION: %d\n", int(dimension_));
  HighsInt num_distances = distances.size();
  if (kDebugReport) printf("num_distances = %d\n", num_distances);
  assert(num_distances == expected_num_distances);
  HighsInt num_col = dimension_ * (dimension_-1);
  std::vector<HighsInt> bi_distance(dimension_ * dimension_);
  i_j_to_k_.clear();
  k_to_i_j_.clear();
  HighsInt k = 0;
  for (HighsInt i = 0; i < dimension_; i++) {
    std::vector<HighsInt> j_to_k(dimension_);
    for (HighsInt j = 0; j < dimension_; j++) {
      if (i == j) {
	j_to_k[j] = -1;
      } else {
	j_to_k[j] = k++;
	k_to_i_j_.push_back(std::make_pair(i, j));
      }
    }
    i_j_to_k_.push_back(j_to_k);
  }
  assert(HighsInt(k_to_i_j_.size()) == num_col);
  assert(HighsInt(i_j_to_k_.size()) == dimension_);
  for (HighsInt i = 0; i < dimension_; i++) 
    assert(HighsInt(i_j_to_k_[i].size()) == dimension_);
  for (HighsInt k = 0; k < num_col; k++) {
    HighsInt iCol = k_to_i_j_[k].first;
    HighsInt jCol = k_to_i_j_[k].second;
    assert(i_j_to_k_[iCol][jCol] == k);
  }

  if (format == "LOWER_DIAG_ROW") {
    k = 0;
    for (HighsInt i = 0; i < dimension_; i++) {
      for (HighsInt j = 0; j < i+1; j++) {
	bi_distance[i+j*dimension_] = distances[k];
	bi_distance[j+i*dimension_] = distances[k];
	k++;
      }
    }
  } else if (format == "FULL_MATRIX") {
    k = 0;
    for (HighsInt i = 0; i < dimension_; i++) {
      for (HighsInt j = 0; j < dimension_; j++) {
	bi_distance[i+j*dimension_] = distances[k];
	k++;
      }
    }
  } 
  lp.num_col_ = num_col;
  lp.num_row_ = 2*dimension_;
  lp.col_cost_.resize(lp.num_col_);
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
  for (HighsInt i = 0; i < dimension_; i++) {
    lp.row_lower_.push_back(1);
    lp.row_upper_.push_back(1);
    for (HighsInt j = 0; j < dimension_; j++) {
      if (i == j) continue;
      HighsInt k = i_j_to_k_[i][j];
      lp.col_cost_[k] = bi_distance[i+j*dimension_];
      lp.a_matrix_.index_.push_back(k);
      lp.a_matrix_.value_.push_back(1);
    }
    lp.a_matrix_.start_.push_back(HighsInt(lp.a_matrix_.index_.size()));
  }
 
  for (HighsInt j = 0; j < dimension_; j++) {
    lp.row_lower_.push_back(1);
    lp.row_upper_.push_back(1);
    for (HighsInt i = 0; i < dimension_; i++) {
      if (i == j) continue;
      HighsInt k = i_j_to_k_[i][j];
      lp.col_cost_[k] = bi_distance[i+j*dimension_];
      lp.a_matrix_.index_.push_back(k);
      lp.a_matrix_.value_.push_back(1);
    }
    lp.a_matrix_.start_.push_back(HighsInt(lp.a_matrix_.index_.size()));
  }
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
}

void TspData::getTours(const double* solution) {
  tours_.clear();
  const HighsInt dimension = dimension_;//(1 + std::sqrt(1+4*num_col)) / 2;
  const HighsInt num_col =  dimension*(dimension-1);
  if (kDebugReport) printf("dimension = %d\n", dimension);
  std::vector<HighsInt> to_node;
  to_node.assign(dimension, -1);
  HighsInt num_nz = 0;
  for (HighsInt k = 0; k < num_col; k++) {
    if (std::abs(solution[k]) < 0.5) continue;
    HighsInt iCol = k_to_i_j_[k].first;
    HighsInt jCol = k_to_i_j_[k].second;
    assert(i_j_to_k_[iCol][jCol] == k);
    assert(to_node[iCol] == -1);
    to_node[iCol] = jCol;
    
    if (kDebugReport) printf("%2d: Edge (%2d, %2d) used\n", int(num_nz), int(iCol), int(jCol));
    num_nz++;
  }
  std::vector<bool> visited;
  visited.assign(dimension, false);
  for (HighsInt k = 0; k < dimension; k++) {
    if (visited[k]) continue;
    HighsInt from_node = k;
    std::vector<HighsInt> tour;
    tour.push_back(from_node);
    visited[from_node] = true;
    for (;;) {
      HighsInt node = to_node[from_node];
      tour.push_back(node);
      if (visited[node]) {
	// Tour identified
	tours_.push_back(tour);
	break;
      } else {
	visited[node] = true;
	from_node = node;
      }
    }
  }
  for (HighsInt k = 0; k < dimension; k++)
    assert(visited[k]);
  if (kDebugReport) {
    HighsInt num_tours = tours_.size();
    for (HighsInt k = 0; k < HighsInt(tours_.size()); k++) {
      printf("Tour %2d: ", int(k));
      std::vector<HighsInt>& tour = tours_[k];
      printf("%2d", int(tour[0]));
      for (HighsInt i=1; i < int(tour.size()); i++)
	printf(" -> %2d", int(tour[i]));
      printf("\n");
    }
  }
  
}

void TspData::addCuts(HighsLp& lp) {
  HighsInt num_tours = tours_.size();
  
  assert(num_tours>1);
  if (num_tours == 0) return;

  for (HighsInt k = 0; k < num_tours; k++) {
    // Add constraint to eliminate this sub-tour
    std::vector<HighsInt>& tour = tours_[k];
    HighsInt tour_num_nodes = tour.size()-1;
    HighsInt from_node = tour[0];
    lp.row_lower_.push_back(-kHighsInf);
    lp.row_upper_.push_back(tour_num_nodes-1);
    for (HighsInt i = 0; i < tour_num_nodes; i++) {
      HighsInt to_node = tour[i+1];
      HighsInt iCol = i_j_to_k_[from_node][to_node];
      lp.a_matrix_.index_.push_back(iCol);
      lp.a_matrix_.value_.push_back(1);
      from_node = to_node;
    }
    lp.a_matrix_.start_.push_back(HighsInt(lp.a_matrix_.index_.size()));
    lp.num_row_++;
  }
}
