#include "HighsRko.h"
#include <chrono>
#include <cstdio>

// ============================================================================
// 1. CONFIG
// ============================================================================

void RKOConfig::print() const {
    printf("\n=== RKO Configuration ===\n");
    printf("Problem: p=%d, alpha=%.2f\n", p, alpha);
    printf("Control: max_time=%.2f, algorithms=%d, runs=%d, ipr=%s\n",
           max_time, num_algorithms, num_runs, enable_ipr ? "on" : "off");
    printf("Elite pool: %d\n", elite_pool_size);
    printf("=========================\n");
}

// ============================================================================
// 2. RANDOM GENERATOR
// ============================================================================

double RandomGenerator::uniform(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng_);
}

int RandomGenerator::uniformInt(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng_);
}

void RkoSolution::randomize() {
    for (size_t i = 0; i < rk.size(); i++) {
        rk[i] = RandomGenerator::instance().uniform(0.0, 1.0);
    }
}

// ============================================================================
// 3. THLP DATA
// ============================================================================

void THLPData::init(int size) {
    n = size;
    W.assign(n, std::vector<double>(n, 0.0));
    C.assign(n, std::vector<double>(n, 0.0));
    O.assign(n, 0.0);
    D.assign(n, 0.0);
    node_names.resize(n);
}

void THLPData::precompute() {
    for (int i = 0; i < n; i++) {
        O[i] = 0.0;
        D[i] = 0.0;
        for (int j = 0; j < n; j++) {
            O[i] += W[i][j];
            D[i] += W[j][i];
        }
    }
}

void THLPData::printSummary() const {
    printf("THLP Data: n=%d, p=%d, alpha=%.2f\n", n, p, alpha);
    double total_flow = 0.0;
    for (std::vector<double>::const_iterator it = O.begin(); it != O.end(); ++it) {
        total_flow += *it;
    }
    printf("  Total flow: %.2f\n", total_flow);
}

void THLPSolution::clear() {
    hubs.clear();
    assignment.clear();
    hub_tree.clear();
    total_cost = 0.0;
    feasible = false;
}

// ============================================================================
// 4. THLP DECODER
// ============================================================================

THLPDecoder::THLPDecoder(const THLPData& data, const RKOConfig& config)
    : data_(data), config_(config) {
    chromosome_size_ = data.n + (data.n - config.p) + (config.p * (config.p - 1) / 2);
}

std::vector<int> THLPDecoder::selectHubs(const std::vector<double>& scores, int p) {
    std::vector<std::pair<double, int> > sorted;
    for (int i = 0; i < (int)scores.size(); i++) {
        sorted.push_back(std::make_pair(scores[i], i));
    }
    std::sort(sorted.begin(), sorted.end(), 
              std::greater<std::pair<double, int> >());
    
    std::vector<int> hubs;
    for (int i = 0; i < p && i < (int)sorted.size(); i++) {
        hubs.push_back(sorted[i].second);
    }
    std::sort(hubs.begin(), hubs.end());
    return hubs;
}

std::vector<int> THLPDecoder::assignNodes(const std::vector<int>& hubs,
                                           const std::vector<double>& keys) {
    int n = data_.n;
    int p = hubs.size();
    std::vector<bool> is_hub(n, false);
    for (std::vector<int>::const_iterator it = hubs.begin(); it != hubs.end(); ++it) {
        is_hub[*it] = true;
    }
    
    std::vector<int> assignment(n, -1);
    for (std::vector<int>::const_iterator it = hubs.begin(); it != hubs.end(); ++it) {
        assignment[*it] = *it;
    }
    
    std::vector<double> hub_ranges(p);
    for (int h = 0; h < p; h++) {
        hub_ranges[h] = (double)(h + 1) / p;
    }
    
    int key_idx = 0;
    for (int node = 0; node < n; node++) {
        if (is_hub[node]) continue;
        if (key_idx >= (int)keys.size()) break;
        
        double key = keys[key_idx++];
        int selected = 0;
        for (int h = 0; h < p; h++) {
            if (key <= hub_ranges[h]) {
                selected = h;
                break;
            }
        }
        assignment[node] = hubs[selected];
    }
    return assignment;
}

std::vector<std::pair<int, int> > THLPDecoder::buildHubTree(
    const std::vector<int>& hubs,
    const std::vector<double>& edge_priorities,
    const std::vector<int>& assignment) {
    
    int p = hubs.size();
    std::vector<std::tuple<double, int, int> > edges;
    
    int edge_idx = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            double priority = (edge_idx < (int)edge_priorities.size()) 
                              ? edge_priorities[edge_idx] 
                              : RandomGenerator::instance().uniform(0.0, 1.0);
            edges.push_back(std::make_tuple(priority, i, j));
            edge_idx++;
        }
    }
    
    std::sort(edges.begin(), edges.end());
    
    std::vector<int> parent(p);
    std::iota(parent.begin(), parent.end(), 0);
    
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    
    std::function<bool(int,int)> unite = [&](int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra != rb) {
            parent[ra] = rb;
            return true;
        }
        return false;
    };
    
    std::vector<std::pair<int, int> > tree;
    int added = 0;
    for (std::vector<std::tuple<double, int, int> >::const_iterator it = edges.begin();
         it != edges.end(); ++it) {
        int u = std::get<1>(*it);
        int v = std::get<2>(*it);
        if (unite(u, v)) {
            tree.push_back(std::make_pair(hubs[u], hubs[v]));
            added++;
            if (added == p - 1) break;
        }
    }
    return tree;
}

double THLPDecoder::calculateCost(const THLPSolution& solution) {
    if (!solution.feasible) return 1e9;
    
    int n = data_.n;
    int p = solution.hubs.size();
    
    std::vector<std::vector<int> > hub_adj(p);
    for (std::vector<std::pair<int, int> >::const_iterator it = solution.hub_tree.begin();
         it != solution.hub_tree.end(); ++it) {
        int u = -1, v = -1;
        for (int i = 0; i < p; i++) {
            if (solution.hubs[i] == it->first) u = i;
            if (solution.hubs[i] == it->second) v = i;
        }
        if (u >= 0 && v >= 0) {
            hub_adj[u].push_back(v);
            hub_adj[v].push_back(u);
        }
    }
    
    std::function<double(int,int)> tree_distance = [&](int start, int end) {
        std::vector<double> dist(p, std::numeric_limits<double>::infinity());
        std::queue<int> q;
        dist[start] = 0.0;
        q.push(start);
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            if (u == end) break;
            for (std::vector<int>::const_iterator it = hub_adj[u].begin();
                 it != hub_adj[u].end(); ++it) {
                int v = *it;
                if (dist[v] > dist[u] + data_.alpha * 
                    data_.C[solution.hubs[u]][solution.hubs[v]]) {
                    dist[v] = dist[u] + data_.alpha * 
                              data_.C[solution.hubs[u]][solution.hubs[v]];
                    q.push(v);
                }
            }
        }
        return dist[end];
    };
    
    double total = 0.0;
    for (int i = 0; i < n; i++) {
        int hub_i = -1;
        for (int h = 0; h < p; h++) {
            if (solution.hubs[h] == solution.assignment[i]) {
                hub_i = h;
                break;
            }
        }
        if (hub_i < 0) continue;
        
        for (int j = 0; j < n; j++) {
            if (data_.W[i][j] == 0) continue;
            
            int hub_j = -1;
            for (int h = 0; h < p; h++) {
                if (solution.hubs[h] == solution.assignment[j]) {
                    hub_j = h;
                    break;
                }
            }
            if (hub_j < 0) continue;
            
            double cost1 = data_.C[i][solution.hubs[hub_i]];
            double cost2 = (hub_i == hub_j) ? 0.0 : tree_distance(hub_i, hub_j);
            double cost3 = data_.C[solution.hubs[hub_j]][j];
            
            total += data_.W[i][j] * (cost1 + cost2 + cost3);
        }
    }
    return total;
}

THLPSolution THLPDecoder::decode(const std::vector<double>& rk,
                                  std::vector<double>& binary_solution) {
    int n = data_.n;
    int p = config_.p;
    
    THLPSolution solution(n);
    
    if (p <= 0 || p > n) {
        solution.feasible = false;
        return solution;
    }
    
    std::vector<double> scores(rk.begin(), rk.begin() + n);
    std::vector<double> assign_keys(rk.begin() + n, rk.begin() + n + (n - p));
    std::vector<double> edge_priorities(rk.begin() + n + (n - p), rk.end());
    
    solution.hubs = selectHubs(scores, p);
    solution.assignment = assignNodes(solution.hubs, assign_keys);
    solution.hub_tree = buildHubTree(solution.hubs, edge_priorities, 
                                      solution.assignment);
    
    // Validate tree
    std::vector<bool> connected(p, false);
    std::vector<std::vector<int> > adj(p);
    for (std::vector<std::pair<int, int> >::const_iterator it = solution.hub_tree.begin();
         it != solution.hub_tree.end(); ++it) {
        int u = -1, v = -1;
        for (int i = 0; i < p; i++) {
            if (solution.hubs[i] == it->first) u = i;
            if (solution.hubs[i] == it->second) v = i;
        }
        if (u >= 0 && v >= 0) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    
    std::queue<int> q;
    connected[0] = true;
    q.push(0);
    int count = 1;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (std::vector<int>::const_iterator it = adj[u].begin();
             it != adj[u].end(); ++it) {
            int v = *it;
            if (!connected[v]) {
                connected[v] = true;
                count++;
                q.push(v);
            }
        }
    }
    solution.feasible = (count == p);
    
    if (solution.feasible) {
        solution.total_cost = calculateCost(solution);
    } else {
        solution.total_cost = 1e9;
    }
    
    binary_solution.assign(n, 0.0);
    for (std::vector<int>::const_iterator it = solution.hubs.begin();
         it != solution.hubs.end(); ++it) {
        binary_solution[*it] = 1.0;
    }
    
    return solution;
}

double THLPDecoder::evaluate(RkoSolution& sol, 
                              std::vector<double>& binary_solution) {
    THLPSolution thlp_sol = decode(sol.rk, binary_solution);
    sol.ofv = thlp_sol.total_cost;
    return sol.ofv;
}

std::function<double(RkoSolution&)> THLPDecoder::createEvaluator(
    std::vector<double>& binary_solution) {
    return std::bind(&THLPDecoder::evaluate, this, std::placeholders::_1, 
                     std::ref(binary_solution));
}

// ============================================================================
// 5. SHAKING - 4 Types 
// ============================================================================

void RKOOptimizer::shakeSolution(RkoSolution& sol, double beta_min, double beta_max) {
    int n = sol.size();
    int intensity = (int)(RandomGenerator::instance().uniform(beta_min, beta_max) * n);
    
    for (int i = 0; i < intensity; i++) {
        int move = RandomGenerator::instance().uniformInt(1, 4);
        int idx = RandomGenerator::instance().uniformInt(0, n - 1);
        
        switch (move) {
            case 1: { // Random
                sol.rk[idx] = RandomGenerator::instance().uniform(0.0, 1.0);
                break;
            }
            case 2: { // Mirror
                sol.rk[idx] = 1.0 - sol.rk[idx];
                break;
            }
            case 3: { // Swap
                int j = RandomGenerator::instance().uniformInt(0, n - 1);
                std::swap(sol.rk[idx], sol.rk[j]);
                break;
            }
            case 4: { // Swap Neighbor
                int j = (idx + 1) % n;
                std::swap(sol.rk[idx], sol.rk[j]);
                break;
            }
        }
    }
}

// ============================================================================
// 6. LOCAL SEARCH VARIANTS
// ============================================================================

void RKOOptimizer::simpleLocalSearch(RkoSolution& sol, 
                                      std::function<double(RkoSolution&)> evaluator,
                                      int iterations) {
    int n = sol.size();
    RkoSolution best = sol;
    best.ofv = evaluator(best);
    
    for (int iter = 0; iter < iterations; iter++) {
        RkoSolution candidate = best;
        int idx = RandomGenerator::instance().uniformInt(0, n - 1);
        candidate.rk[idx] += RandomGenerator::instance().uniform(-config_.ls_step_size, 
                                                                  config_.ls_step_size);
        candidate.rk[idx] = std::max(0.0, std::min(1.0, candidate.rk[idx]));
        candidate.ofv = evaluator(candidate);
        
        if (candidate.ofv < best.ofv) {
            best = candidate;
        }
    }
    sol = best;
}

void RKOOptimizer::strongLocalSearch(RkoSolution& sol, 
                                      std::function<double(RkoSolution&)> evaluator,
                                      int iterations) {
    int n = sol.size();
    RkoSolution best = sol;
    best.ofv = evaluator(best);
    
    for (int iter = 0; iter < iterations; iter++) {
        RkoSolution candidate = best;
        int num_changes = RandomGenerator::instance().uniformInt(1, 3);
        for (int k = 0; k < num_changes; k++) {
            int idx = RandomGenerator::instance().uniformInt(0, n - 1);
            candidate.rk[idx] += RandomGenerator::instance().uniform(-0.1, 0.1);
            candidate.rk[idx] = std::max(0.0, std::min(1.0, candidate.rk[idx]));
        }
        candidate.ofv = evaluator(candidate);
        
        if (candidate.ofv < best.ofv) {
            best = candidate;
        }
    }
    sol = best;
}

void RKOOptimizer::swapLocalSearch(RkoSolution& sol,
                                    std::function<double(RkoSolution&)> evaluator) {
    int n = sol.size();
    RkoSolution best = sol;
    best.ofv = evaluator(best);
    
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), RandomGenerator::instance().getEngine());
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            RkoSolution candidate = best;
            std::swap(candidate.rk[order[i]], candidate.rk[order[j]]);
            candidate.ofv = evaluator(candidate);
            
            if (candidate.ofv < best.ofv) {
                best = candidate;
            }
        }
    }
    sol = best;
}

void RKOOptimizer::mirrorLocalSearch(RkoSolution& sol,
                                      std::function<double(RkoSolution&)> evaluator) {
    int n = sol.size();
    RkoSolution best = sol;
    best.ofv = evaluator(best);
    
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), RandomGenerator::instance().getEngine());
    
    for (int i = 0; i < n; i++) {
        RkoSolution candidate = best;
        int idx = order[i];
        candidate.rk[idx] = 1.0 - candidate.rk[idx];
        candidate.ofv = evaluator(candidate);
        
        if (candidate.ofv < best.ofv) {
            best = candidate;
        }
    }
    sol = best;
}

// ============================================================================
// 7. FAREY LOCAL SEARCH
// ============================================================================

static void fareySequence(int num, std::vector<double>& F) {
    int a = 0, b = 1, c = 1, d = num;
    F.push_back((double)a / b);
    while (c <= num) {
        int k = (num + b) / d;
        int temp_a = a, temp_b = b;
        a = c;
        b = d;
        c = k * c - temp_a;
        d = k * d - temp_b;
        F.push_back((double)a / b);
    }
}

void RKOOptimizer::fareyLocalSearch(RkoSolution& sol,
                                     std::function<double(RkoSolution&)> evaluator) {
    std::vector<double> F;
    fareySequence(7, F);
    
    int n = sol.size();
    RkoSolution best = sol;
    best.ofv = evaluator(best);
    
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), RandomGenerator::instance().getEngine());
    
    for (int i = 0; i < n; i++) {
        int idx = order[i];
        double best_val = best.rk[idx];
        double best_ofv = best.ofv;
        
        for (size_t j = 0; j < F.size() - 1; j++) {
            RkoSolution candidate = best;
            candidate.rk[idx] = RandomGenerator::instance().uniform(F[j], F[j + 1]);
            candidate.ofv = evaluator(candidate);
            
            if (candidate.ofv < best_ofv) {
                best_ofv = candidate.ofv;
                best_val = candidate.rk[idx];
            }
        }
        best.rk[idx] = best_val;
        best.ofv = best_ofv;
    }
    sol = best;
}

// ============================================================================
// 8. NELDER-MEAD
// ============================================================================

RkoSolution RKOOptimizer::blend(const RkoSolution& a, const RkoSolution& b, double rho) {
    int n = a.size();
    RkoSolution result(n);
    for (int i = 0; i < n; i++) {
        result.rk[i] = a.rk[i] + rho * (b.rk[i] - a.rk[i]);
        result.rk[i] = std::max(0.0, std::min(1.0, result.rk[i]));
    }
    return result;
}

void RKOOptimizer::nelderMeadSearch(RkoSolution& sol,
                                     std::function<double(RkoSolution&)> evaluator,
                                     const std::vector<RkoSolution>& pool) {
    int n = sol.size();
    
    if (pool.size() < 2) {
        strongLocalSearch(sol, evaluator, config_.ls_strong_iterations);
        return;
    }
    
    RkoSolution x1 = sol;
    x1.ofv = evaluator(x1);
    
    int idx2 = RandomGenerator::instance().uniformInt(0, pool.size() - 1);
    int idx3;
    do {
        idx3 = RandomGenerator::instance().uniformInt(0, pool.size() - 1);
    } while (idx2 == idx3);
    
    RkoSolution x2 = pool[idx2];
    RkoSolution x3 = pool[idx3];
    
    std::vector<RkoSolution> simplex;
    simplex.push_back(x1);
    simplex.push_back(x2);
    simplex.push_back(x3);
    std::sort(simplex.begin(), simplex.end(),
              [](const RkoSolution& a, const RkoSolution& b) {
                  return a.ofv < b.ofv;
              });
    
    int max_iter = (int)(n * 0.135);
    max_iter = std::max(3, max_iter);
    
    for (int iter = 0; iter < max_iter; iter++) {
        RkoSolution best = simplex[0];
        RkoSolution middle = simplex[1];
        RkoSolution worst = simplex[2];
        
        RkoSolution centroid = blend(best, middle, 0.5);
        centroid.ofv = evaluator(centroid);
        
        RkoSolution reflected = blend(centroid, worst, -1.0);
        reflected.ofv = evaluator(reflected);
        
        if (reflected.ofv < best.ofv) {
            RkoSolution expanded = blend(reflected, centroid, -1.0);
            expanded.ofv = evaluator(expanded);
            simplex[2] = (expanded.ofv < reflected.ofv) ? expanded : reflected;
        } else if (reflected.ofv < middle.ofv) {
            simplex[2] = reflected;
        } else {
            RkoSolution contracted;
            if (reflected.ofv < worst.ofv) {
                contracted = blend(reflected, centroid, 1.0);
            } else {
                contracted = blend(centroid, worst, 1.0);
            }
            contracted.ofv = evaluator(contracted);
            
            if (contracted.ofv < worst.ofv) {
                simplex[2] = contracted;
            } else {
                for (int i = 1; i < 3; i++) {
                    simplex[i] = blend(best, simplex[i], 1.0);
                    simplex[i].ofv = evaluator(simplex[i]);
                }
            }
        }
        
        std::sort(simplex.begin(), simplex.end(),
                  [](const RkoSolution& a, const RkoSolution& b) {
                      return a.ofv < b.ofv;
                  });
    }
    
    sol = simplex[0];
    sol.ofv = evaluator(sol);
}

// ============================================================================
// 9. RVND - Randomized Variable Neighborhood Descent 
// ============================================================================

void RKOOptimizer::rvnd(RkoSolution& sol, std::function<double(RkoSolution&)> evaluator) {
    // 4 neighborhoods: Swap, Farey, Mirror, Simple
    std::vector<int> neighbors = {1, 2, 3, 4};
    
    while (!neighbors.empty()) {
        int idx = RandomGenerator::instance().uniformInt(0, neighbors.size() - 1);
        int nbr = neighbors[idx];
        
        RkoSolution candidate = sol;
        
        switch (nbr) {
            case 1: swapLocalSearch(candidate, evaluator); break;
            case 2: fareyLocalSearch(candidate, evaluator); break;
            case 3: mirrorLocalSearch(candidate, evaluator); break;
            case 4: simpleLocalSearch(candidate, evaluator, 30); break;
        }
        
        if (candidate.ofv < sol.ofv) {
            sol = candidate;
            // Restart neighborhood list
            neighbors.clear();
            neighbors.push_back(1);
            neighbors.push_back(2);
            neighbors.push_back(3);
            neighbors.push_back(4);
        } else {
            neighbors.erase(neighbors.begin() + idx);
        }
    }
}

// ============================================================================
// 10. IPR
// ============================================================================

RkoSolution RKOOptimizer::iprPath(const RkoSolution& atual, const RkoSolution& guia,
                                   std::function<double(RkoSolution&)> evaluator) {
    int n = atual.size();
    int blockSize = std::max(1, (int)(n * config_.ipr_block_ratio));
    int numBlock = std::max(1, n / blockSize);
    
    std::vector<int> fixedBlock(numBlock, 1);
    int dist = 0;
    
    for (int i = 0; i < n; i++) {
        if (atual.rk[i] != guia.rk[i]) {
            fixedBlock[i % numBlock] = 0;
            dist++;
        }
    }
    
    if (dist == 0) return atual;
    
    RkoSolution sCurrent = atual;
    RkoSolution bestPath = atual;
    bestPath.ofv = evaluator(bestPath);
    int direction = config_.ipr_direction;
    int maxIter = numBlock - 1;
    int iter = 0;
    
    while (iter < maxIter) {
        iter++;
        int bestBlock = -1;
        RkoSolution bestIter = sCurrent;
        bestIter.ofv = std::numeric_limits<double>::infinity();
        
        for (int i = 0; i < numBlock; i++) {
            if (fixedBlock[i] == 0) {
                RkoSolution sViz = sCurrent;
                int start = i * blockSize;
                int end = std::min(start + blockSize, n);
                
                for (int k = start; k < end; k++) {
                    if (direction == 1) {
                        sViz.rk[k] = guia.rk[k];
                    } else {
                        sViz.rk[k] = std::max(0.0, std::min(1.0 - guia.rk[k], 0.999999));
                    }
                }
                sViz.ofv = evaluator(sViz);
                
                if (sViz.ofv < bestIter.ofv) {
                    bestIter = sViz;
                    bestBlock = i;
                }
            }
        }
        
        if (bestIter.ofv < bestPath.ofv) {
            bestPath = bestIter;
            if (config_.verbose) {
                printf("      IPR improved: %.2f\n", bestPath.ofv);
            }
        }
        
        if (bestBlock >= 0) {
            sCurrent = bestIter;
            fixedBlock[bestBlock] = 1;
        }
    }
    
    return bestPath;
}

void RKOOptimizer::runIPR(std::function<double(RkoSolution&)> evaluator) {
    if (!config_.enable_ipr || elite_pool_.size() < 2) return;
    
    int idx1 = RandomGenerator::instance().uniformInt(0, elite_pool_.size() - 1);
    int idx2;
    do {
        idx2 = RandomGenerator::instance().uniformInt(0, elite_pool_.size() - 1);
    } while (idx1 == idx2);
    
    RkoSolution result = iprPath(elite_pool_[idx1], elite_pool_[idx2], evaluator);
    updateElitePool(result);
}

// ============================================================================
// 11. ELITE POOL
// ============================================================================

void RKOOptimizer::updateElitePool(const RkoSolution& sol) {
    // Check similarity
    for (std::vector<RkoSolution>::iterator it = elite_pool_.begin(); 
         it != elite_pool_.end(); ++it) {
        double diff = 0.0;
        for (size_t i = 0; i < sol.rk.size() && i < it->rk.size(); i++) {
            diff += std::abs(sol.rk[i] - it->rk[i]);
        }
        diff /= sol.rk.size();
        if (diff < 0.01) return;
    }
    
    elite_pool_.push_back(sol);
    std::sort(elite_pool_.begin(), elite_pool_.end(),
              [](const RkoSolution& a, const RkoSolution& b) {
                  return a.ofv < b.ofv;
              });
    
    if ((int)elite_pool_.size() > config_.elite_pool_size) {
        elite_pool_.resize(config_.elite_pool_size);
    }
}

// ============================================================================
// 12. SA ALGORITHM
// ============================================================================

bool SAAlgorithm::run(RkoSolution& sol,
                       std::function<double(RkoSolution&)> evaluator,
                       const RKOConfig& config,
                       bool verbose) {
    double T = config.sa_T0;
    int iter = 0;
    int improvements = 0;
    int accepted = 0;
    
    sol.ofv = evaluator(sol);
    RkoSolution best = sol;
    RkoSolution current = sol;
    
    if (verbose) {
        printf("    SA: T0=%.0f, alpha=%.3f, max_iter=%d\n",
               config.sa_T0, config.sa_alpha, config.sa_max_iter);
        printf("    beta=[%.2f,%.2f]\n", config.sa_beta_min, config.sa_beta_max);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    
    while (elapsed < config.max_time) {
        for (int i = 0; i < config.sa_max_iter && elapsed < config.max_time; i++) {
            iter++;
            RkoSolution neighbor = current;
            
            // Use RKOOptimizer's shake method
            RKOOptimizer tmp;
            tmp.shakeSolution(neighbor, config.sa_beta_min, config.sa_beta_max);
            neighbor.ofv = evaluator(neighbor);
            
            double delta = neighbor.ofv - current.ofv;
            if (delta < 0 || RandomGenerator::instance().uniform(0.0, 1.0) < 
                std::exp(-delta / T)) {
                current = neighbor;
                accepted++;
                if (current.ofv < best.ofv) {
                    best = current;
                    improvements++;
                }
            }
            
            elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
        }
        T *= config.sa_alpha;
        
        // Progress every ~20%
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | T=%.0f | best=%.2f\n",
                   elapsed / config.max_time * 100, T, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    SA: %d iters, %d improvements, %d accepted, best=%.2f\n",
               iter, improvements, accepted, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 13. ILS ALGORITHM
// ============================================================================

bool ILSAlgorithm::run(RkoSolution& sol,
                        std::function<double(RkoSolution&)> evaluator,
                        const RKOConfig& config,
                        bool verbose) {
    sol.ofv = evaluator(sol);
    RkoSolution best = sol;
    int improvements = 0;
    int iter = 0;
    
    if (verbose) {
        printf("    ILS: beta=[%.2f,%.2f]\n",
               config.ils_beta_min, config.ils_beta_max);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    
    while (elapsed < config.max_time) {
        iter++;
        RkoSolution candidate = best;
        
        RKOOptimizer tmp;
        tmp.shakeSolution(candidate, config.ils_beta_min, config.ils_beta_max);
        candidate.ofv = evaluator(candidate);
        
        // Local search (RVND)
        tmp.rvnd(candidate, evaluator);
        
        if (candidate.ofv < best.ofv) {
            best = candidate;
            improvements++;
            if (verbose && config.verbose) {
                printf("      ILS improved: %.2f (iter %d)\n", best.ofv, iter);
            }
        }
        
        elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    ILS: %d iters, %d improvements, best=%.2f\n",
               iter, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 14. VNS ALGORITHM
// ============================================================================

bool VNSAlgorithm::run(RkoSolution& sol,
                        std::function<double(RkoSolution&)> evaluator,
                        const RKOConfig& config,
                        bool verbose) {
    sol.ofv = evaluator(sol);
    RkoSolution best = sol;
    RkoSolution current = sol;
    int improvements = 0;
    int iter = 0;
    
    if (verbose) {
        printf("    VNS: k_max=%d, beta_min=%.3f\n",
               config.vns_k_max, config.vns_beta_min);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        int k = 1;
        while (k <= config.vns_k_max && elapsed < config.max_time) {
            iter++;
            RkoSolution candidate = current;
            
            double beta = RandomGenerator::instance().uniform(
                k * config.vns_beta_min, (k + 1) * config.vns_beta_min);
            tmp.shakeSolution(candidate, beta, beta);
            candidate.ofv = evaluator(candidate);
            tmp.rvnd(candidate, evaluator);
            
            if (candidate.ofv < current.ofv) {
                current = candidate;
                if (current.ofv < best.ofv) {
                    best = current;
                    improvements++;
                }
                k = 1;
            } else {
                k++;
            }
            
            elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
        }
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    VNS: %d iters, %d improvements, best=%.2f\n",
               iter, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 15. LNS ALGORITHM (with Farey)
// ============================================================================

bool LNSAlgorithm::run(RkoSolution& sol,
                        std::function<double(RkoSolution&)> evaluator,
                        const RKOConfig& config,
                        bool verbose) {
    std::vector<double> F;
    fareySequence(config.lns_farey_level, F);
    
    int n = sol.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    
    sol.ofv = evaluator(sol);
    RkoSolution best = sol;
    RkoSolution current = sol;
    
    if (verbose) {
        printf("    LNS: beta=[%.3f,%.3f], farey_level=%d\n",
               config.lns_beta_min, config.lns_beta_max, config.lns_farey_level);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int iter = 0;
    double T = config.lns_T0;
    int improvements = 0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        iter++;
        RkoSolution candidate = current;
        int intensity = RandomGenerator::instance().uniformInt(
            (int)(config.lns_beta_min * n), (int)(config.lns_beta_max * n));
        std::shuffle(order.begin(), order.end(), 
                     RandomGenerator::instance().getEngine());
        
        for (int k = 0; k < intensity && k < n; k++) {
            int pos = order[k];
            double best_val = current.rk[pos];
            double best_ofv = std::numeric_limits<double>::infinity();
            
            for (size_t j = 0; j < F.size() - 1; j++) {
                double val = RandomGenerator::instance().uniform(F[j], F[j + 1]);
                candidate.rk[pos] = val;
                double ofv = evaluator(candidate);
                if (ofv < best_ofv) {
                    best_ofv = ofv;
                    best_val = val;
                }
            }
            candidate.rk[pos] = best_val;
            candidate.ofv = best_ofv;
        }
        
        tmp.rvnd(candidate, evaluator);
        
        double delta = candidate.ofv - current.ofv;
        if (delta <= 0 || RandomGenerator::instance().uniform(0.0, 1.0) < 
            std::exp(-delta / T)) {
            current = candidate;
            if (current.ofv < best.ofv) {
                best = current;
                improvements++;
            }
        }
        
        T *= config.lns_alpha;
        elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    LNS: %d iters, %d improvements, best=%.2f\n",
               iter, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 16. GRASP ALGORITHM
// ============================================================================

static void constructiveGRASP(RkoSolution& sol,
                               std::function<double(RkoSolution&)> evaluator,
                               double h) {
    int n = sol.size();
    std::vector<int> unfixed(n);
    std::iota(unfixed.begin(), unfixed.end(), 0);
    
    double intensity = RandomGenerator::instance().uniform(0.3, 0.7);
    int iterations = (int)(n * intensity);
    
    for (int it = 0; it < iterations && !unfixed.empty(); it++) {
        int kMax = std::max(2, (int)(unfixed.size() * 0.1));
        kMax = std::min(kMax, (int)unfixed.size());
        std::shuffle(unfixed.begin(), unfixed.end(), 
                     RandomGenerator::instance().getEngine());
        
        std::vector<double> best_val(n, 0.0);
        std::vector<double> best_ofv(n, std::numeric_limits<double>::infinity());
        
        for (int k = 0; k < kMax; k++) {
            int idx = unfixed[k];
            RkoSolution aux = sol;
            double best_z = aux.rk[idx];
            double best_f = evaluator(aux);
            
            std::vector<double> values;
            values.push_back(aux.rk[idx]);
            for (int j = 0; j < (int)(1.0 / h) + 1; j += 2) {
                if (aux.rk[idx] + j * h >= 0 && aux.rk[idx] + j * h < 1)
                    values.push_back(aux.rk[idx] + j * h);
                if (aux.rk[idx] - j * h >= 0 && aux.rk[idx] - j * h < 1)
                    values.push_back(aux.rk[idx] - j * h);
            }
            
            std::shuffle(values.begin(), values.end(),
                         RandomGenerator::instance().getEngine());
            int q = std::min((int)values.size(), 
                             (int)std::ceil(std::log2(1.0 / h)) + 1);
            
            for (int j = 0; j < q; j++) {
                aux.rk[idx] = values[j];
                double f = evaluator(aux);
                if (f < best_f) {
                    best_f = f;
                    best_z = values[j];
                }
            }
            
            best_val[idx] = best_z;
            best_ofv[idx] = best_f;
        }
        
        int best_idx = -1;
        double min_f = std::numeric_limits<double>::infinity();
        for (int k = 0; k < kMax; k++) {
            int idx = unfixed[k];
            if (best_ofv[idx] < min_f) {
                min_f = best_ofv[idx];
                best_idx = idx;
            }
        }
        
        if (best_idx >= 0) {
            sol.rk[best_idx] = best_val[best_idx];
            sol.ofv = best_ofv[best_idx];
            unfixed.erase(std::remove(unfixed.begin(), unfixed.end(), best_idx),
                          unfixed.end());
        }
    }
    sol.ofv = evaluator(sol);
}

bool GRASPAlgorithm::run(RkoSolution& sol,
                          std::function<double(RkoSolution&)> evaluator,
                          const RKOConfig& config,
                          bool verbose) {
    sol.ofv = evaluator(sol);
    RkoSolution best = sol;
    RkoSolution current = sol;
    int improvements = 0;
    int iter = 0;
    
    if (verbose) {
        printf("    GRASP: hs=%.4f, he=%.6f\n",
               config.grasp_hs, config.grasp_he);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        iter++;
        double h = config.grasp_hs;
        while (h >= config.grasp_he && elapsed < config.max_time) {
            RkoSolution candidate = current;
            constructiveGRASP(candidate, evaluator, h);
            tmp.rvnd(candidate, evaluator);
            
            if (candidate.ofv < best.ofv) {
                best = candidate;
                improvements++;
            }
            
            if (candidate.ofv < current.ofv) {
                current = candidate;
            } else if (RandomGenerator::instance().uniform(0.0, 1.0) < 
                       std::exp(-(candidate.ofv - current.ofv) / 
                                (100 - 100 * elapsed / config.max_time))) {
                current = candidate;
            }
            
            h /= 2;
            elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
        }
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    GRASP: %d iters, %d improvements, best=%.2f\n",
               iter, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 17. GA ALGORITHM
// ============================================================================

bool GAAlgorithm::run(RkoSolution& sol,
                       std::function<double(RkoSolution&)> evaluator,
                       const RKOConfig& config,
                       bool verbose) {
    int n = sol.size();
    int pop_size = config.ga_population_size;
    
    std::vector<RkoSolution> population(pop_size, RkoSolution(n));
    for (int i = 0; i < pop_size; i++) {
        population[i].randomize();
        population[i].ofv = evaluator(population[i]);
    }
    std::sort(population.begin(), population.end(),
              [](const RkoSolution& a, const RkoSolution& b) {
                  return a.ofv < b.ofv;
              });
    
    RkoSolution best = population[0];
    int improvements = 0;
    int gen = 0;
    
    if (verbose) {
        printf("    GA: pop=%d, init=%.2f\n", pop_size, best.ofv);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        gen++;
        std::vector<RkoSolution> offspring(pop_size, RkoSolution(n));
        
        for (int i = 0; i < pop_size; i += 2) {
            std::function<RkoSolution()> tournament = [&]() {
                int best_idx = RandomGenerator::instance().uniformInt(0, pop_size - 1);
                for (int j = 1; j < config.ga_tournament_size; j++) {
                    int idx = RandomGenerator::instance().uniformInt(0, pop_size - 1);
                    if (population[idx].ofv < population[best_idx].ofv) {
                        best_idx = idx;
                    }
                }
                return population[best_idx];
            };
            
            RkoSolution p1 = tournament();
            RkoSolution p2 = tournament();
            
            offspring[i] = p1;
            if (i + 1 < pop_size) offspring[i + 1] = p2;
            
            if (RandomGenerator::instance().uniform(0.0, 1.0) < config.ga_crossover_prob) {
                for (int j = 0; j < n; j++) {
                    if (RandomGenerator::instance().uniform(0.0, 1.0) < 0.5) {
                        offspring[i].rk[j] = p2.rk[j];
                        if (i + 1 < pop_size) {
                            offspring[i + 1].rk[j] = p1.rk[j];
                        }
                    }
                    if (RandomGenerator::instance().uniform(0.0, 1.0) < config.ga_mutation_prob) {
                        offspring[i].rk[j] = RandomGenerator::instance().uniform(0.0, 1.0);
                    }
                    if (i + 1 < pop_size && 
                        RandomGenerator::instance().uniform(0.0, 1.0) < config.ga_mutation_prob) {
                        offspring[i + 1].rk[j] = RandomGenerator::instance().uniform(0.0, 1.0);
                    }
                }
            }
            
            offspring[i].ofv = evaluator(offspring[i]);
            if (offspring[i].ofv < best.ofv) {
                best = offspring[i];
                improvements++;
            }
            
            if (i + 1 < pop_size) {
                offspring[i + 1].ofv = evaluator(offspring[i + 1]);
                if (offspring[i + 1].ofv < best.ofv) {
                    best = offspring[i + 1];
                    improvements++;
                }
            }
        }
        
        // Elitism
        std::sort(offspring.begin(), offspring.end(),
                  [](const RkoSolution& a, const RkoSolution& b) {
                      return a.ofv < b.ofv;
                  });
        offspring[0] = population[0];
        population = offspring;
        elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    GA: %d gens, %d improvements, best=%.2f\n",
               gen, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = gen;
    return true;
}

// ============================================================================
// 18. PSO ALGORITHM
// ============================================================================

bool PSOAlgorithm::run(RkoSolution& sol,
                        std::function<double(RkoSolution&)> evaluator,
                        const RKOConfig& config,
                        bool verbose) {
    int n = sol.size();
    int size = config.pso_particle_size;
    
    std::vector<RkoSolution> X(size, RkoSolution(n));
    std::vector<RkoSolution> Pbest(size, RkoSolution(n));
    std::vector<std::vector<double> > V(size, std::vector<double>(n, 0.0));
    RkoSolution Gbest(n);
    Gbest.ofv = std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < size; i++) {
        X[i].randomize();
        X[i].ofv = evaluator(X[i]);
        Pbest[i] = X[i];
        for (int j = 0; j < n; j++) {
            V[i][j] = RandomGenerator::instance().uniform(-0.1, 0.1);
        }
        if (X[i].ofv < Gbest.ofv) Gbest = X[i];
    }
    
    if (verbose) {
        printf("    PSO: particles=%d, init=%.2f\n", size, Gbest.ofv);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int iter = 0;
    int improvements = 0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        iter++;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < n; j++) {
                double r1 = RandomGenerator::instance().uniform(0.0, 1.0);
                double r2 = RandomGenerator::instance().uniform(0.0, 1.0);
                V[i][j] = config.pso_w * V[i][j] +
                          config.pso_c1 * r1 * (Pbest[i].rk[j] - X[i].rk[j]) +
                          config.pso_c2 * r2 * (Gbest.rk[j] - X[i].rk[j]);
                
                double new_val = X[i].rk[j] + V[i][j];
                if (new_val >= 0.0 && new_val <= 1.0) {
                    X[i].rk[j] = new_val;
                } else {
                    V[i][j] = 0;
                }
            }
            
            X[i].ofv = evaluator(X[i]);
            if (X[i].ofv < Pbest[i].ofv) {
                Pbest[i] = X[i];
            }
            if (X[i].ofv < Gbest.ofv) {
                Gbest = X[i];
                improvements++;
            }
        }
        
        // Local search on random particle
        int idx = RandomGenerator::instance().uniformInt(0, size - 1);
        tmp.rvnd(Pbest[idx], evaluator);
        if (Pbest[idx].ofv < Gbest.ofv) {
            Gbest = Pbest[idx];
            improvements++;
        }
        
        elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, Gbest.ofv);
        }
    }
    
    if (verbose) {
        printf("    PSO: %d iters, %d improvements, best=%.2f\n",
               iter, improvements, Gbest.ofv);
    }
    
    sol = Gbest;
    sol.ofv = evaluator(sol);
    sol.generation = iter;
    return true;
}

// ============================================================================
// 19. BRKGA ALGORITHM
// ============================================================================

bool BRKGAAlgorithm::run(RkoSolution& sol,
                          std::function<double(RkoSolution&)> evaluator,
                          const RKOConfig& config,
                          bool verbose) {
    int n = sol.size();
    int pop_size = config.brkga_population_size;
    int elite_size = (int)(pop_size * config.brkga_pe);
    int mutant_size = (int)(pop_size * config.brkga_pm);
    
    std::vector<RkoSolution> population(pop_size, RkoSolution(n));
    for (int i = 0; i < pop_size; i++) {
        population[i].randomize();
        population[i].ofv = evaluator(population[i]);
    }
    std::sort(population.begin(), population.end(),
              [](const RkoSolution& a, const RkoSolution& b) {
                  return a.ofv < b.ofv;
              });
    
    RkoSolution best = population[0];
    int improvements = 0;
    int gen = 0;
    
    if (verbose) {
        printf("    BRKGA: pop=%d, init=%.2f\n", pop_size, best.ofv);
    }
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed = 0.0;
    int last_print = 0;
    RKOOptimizer tmp;
    
    while (elapsed < config.max_time) {
        gen++;
        std::vector<RkoSolution> new_pop(pop_size, RkoSolution(n));
        
        // Elite
        for (int i = 0; i < elite_size; i++) {
            new_pop[i] = population[i];
        }
        
        // Crossover
        for (int i = elite_size; i < pop_size - mutant_size; i++) {
            int elite_parent = RandomGenerator::instance().uniformInt(0, elite_size - 1);
            int non_elite_parent = RandomGenerator::instance().uniformInt(elite_size, pop_size - 1);
            
            for (int j = 0; j < n; j++) {
                if (RandomGenerator::instance().uniform(0.0, 1.0) < config.brkga_rhoe) {
                    new_pop[i].rk[j] = population[elite_parent].rk[j];
                } else {
                    new_pop[i].rk[j] = population[non_elite_parent].rk[j];
                }
            }
            new_pop[i].ofv = evaluator(new_pop[i]);
        }
        
        // Mutants
        for (int i = pop_size - mutant_size; i < pop_size; i++) {
            new_pop[i].randomize();
            new_pop[i].ofv = evaluator(new_pop[i]);
        }
        
        // Local search on elites
        for (int i = 0; i < std::min(2, elite_size); i++) {
            tmp.rvnd(new_pop[i], evaluator);
        }
        
        population = new_pop;
        std::sort(population.begin(), population.end(),
                  [](const RkoSolution& a, const RkoSolution& b) {
                      return a.ofv < b.ofv;
                  });
        
        if (population[0].ofv < best.ofv) {
            best = population[0];
            improvements++;
        }
        
        elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        
        if (verbose && config.verbose && (int)(elapsed / config.max_time * 5) > last_print) {
            last_print = (int)(elapsed / config.max_time * 5);
            printf("      progress %.0f%% | best=%.2f\n",
                   elapsed / config.max_time * 100, best.ofv);
        }
    }
    
    if (verbose) {
        printf("    BRKGA: %d gens, %d improvements, best=%.2f\n",
               gen, improvements, best.ofv);
    }
    
    sol = best;
    sol.ofv = evaluator(sol);
    sol.generation = gen;
    return true;
}

// ============================================================================
// 20. RKO OPTIMIZER
// ============================================================================

RKOOptimizer::RKOOptimizer() : decoder_(NULL) {}

RKOOptimizer::RKOOptimizer(const RKOConfig& config) : config_(config), decoder_(NULL) {}

RKOOptimizer::~RKOOptimizer() {
    if (decoder_ != NULL) {
        delete decoder_;
        decoder_ = NULL;
    }
}

RkoSolution RKOOptimizer::createInitialSolution(int chromosome_size) {
    RkoSolution sol(chromosome_size);
    sol.randomize();
    return sol;
}

bool RKOOptimizer::runSingleAlgorithm(RKOAlgorithm& algorithm,
                                       RkoSolution& sol,
                                       std::function<double(RkoSolution&)> evaluator) {
    sol.ofv = evaluator(sol);
    bool success = algorithm.run(sol, evaluator, config_, config_.verbose);
    sol.ofv = evaluator(sol);
    return success;
}

void RKOOptimizer::updateBest(const RkoSolution& sol) {
    if (sol.ofv < best_rk_solution_.ofv) {
        best_rk_solution_ = sol;
        if (decoder_ != NULL) {
            best_solution_ = decoder_->decode(sol.rk, best_binary_solution_);
        }
    }
}

bool RKOOptimizer::solveTHLP(const THLPData& data, std::vector<double>& solution) {
    printf("\n=== RKO Optimizer ===\n");
    data.printSummary();
    
    config_.p = data.p;
    config_.alpha = data.alpha;
    
    if (decoder_ != NULL) {
        delete decoder_;
        decoder_ = NULL;
    }
    decoder_ = new THLPDecoder(data, config_);
    
    int chromosome_size = decoder_->getChromosomeSize();
    best_rk_solution_.ofv = std::numeric_limits<double>::infinity();
    elite_pool_.clear();
    
    printf("Chromosome: %d, Algorithms: %d, Runs: %d, IPR: %s\n\n",
           chromosome_size, config_.num_algorithms, config_.num_runs,
           config_.enable_ipr ? "on" : "off");
    
    // Algorithm list
    std::vector<std::pair<RKOAlgorithm*, std::string> > algorithms;
    algorithms.push_back(std::make_pair(&sa_, "SA"));
    algorithms.push_back(std::make_pair(&ils_, "ILS"));
    algorithms.push_back(std::make_pair(&vns_, "VNS"));
    algorithms.push_back(std::make_pair(&lns_, "LNS"));
    algorithms.push_back(std::make_pair(&grasp_, "GRASP"));
    algorithms.push_back(std::make_pair(&ga_, "GA"));
    algorithms.push_back(std::make_pair(&pso_, "PSO"));
    algorithms.push_back(std::make_pair(&brkga_, "BRKGA"));
    
    std::vector<std::pair<std::string, double> > results;
    int algo_count = std::min(config_.num_algorithms, (int)algorithms.size());
    
    for (int a = 0; a < algo_count; a++) {
        std::string name = algorithms[a].second;
        printf("[%d/%d] %s\n", a + 1, algo_count, name.c_str());
        
        double best_for_algo = std::numeric_limits<double>::infinity();
        double total_time = 0.0;
        int feasible_runs = 0;
        
        for (int r = 0; r < config_.num_runs; r++) {
            RkoSolution sol = createInitialSolution(chromosome_size);
            std::vector<double> binary_sol;
            std::function<double(RkoSolution&)> evaluator = decoder_->createEvaluator(binary_sol);
            sol.ofv = evaluator(sol);
            
            printf("  Run %d: initial=%.2f", r + 1, sol.ofv);
            
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            runSingleAlgorithm(*algorithms[a].first, sol, evaluator);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(end - start).count();
            total_time += elapsed;
            
            // Check feasible
            std::vector<double> temp_bin;
            THLPSolution check = decoder_->decode(sol.rk, temp_bin);
            if (check.feasible) {
                feasible_runs++;
                printf(" → %.2f (%.2fs)\n", sol.ofv, elapsed);
            } else {
                printf(" → infeasible (%.2fs)\n", elapsed);
            }
            
            if (sol.ofv < best_for_algo) {
                best_for_algo = sol.ofv;
            }
            
            if (sol.ofv < best_rk_solution_.ofv) {
                best_rk_solution_ = sol;
                best_binary_solution_ = binary_sol;
                best_solution_ = decoder_->decode(sol.rk, best_binary_solution_);
            }
            
            updateElitePool(sol);
        }
        
        results.push_back(std::make_pair(name, best_for_algo));
        printf("  %s best: %.2f (avg %.2fs, %d/%d feasible)\n\n",
               name.c_str(), best_for_algo, total_time / config_.num_runs,
               feasible_runs, config_.num_runs);
    }
    
    // Run IPR
    if (config_.enable_ipr && elite_pool_.size() >= 2) {
        printf("Running IPR on elite pool (%zu solutions)...\n", elite_pool_.size());
        std::vector<double> dummy;
        std::function<double(RkoSolution&)> evaluator = decoder_->createEvaluator(dummy);
        runIPR(evaluator);
        
        if (!elite_pool_.empty() && elite_pool_[0].ofv < best_rk_solution_.ofv) {
            best_rk_solution_ = elite_pool_[0];
            best_solution_ = decoder_->decode(best_rk_solution_.rk, best_binary_solution_);
        }
        printf("  Elite pool best: %.2f\n", elite_pool_[0].ofv);
    }
    
    // Sort results
    std::sort(results.begin(), results.end(),
              [](const std::pair<std::string, double>& a, 
                 const std::pair<std::string, double>& b) { 
                  return a.second < b.second; 
              });
    
    printf("\n=== Results ===\n");
    for (size_t i = 0; i < results.size(); i++) {
        printf("%2lu. %-8s: %.2f\n", i + 1, 
               results[i].first.c_str(), results[i].second);
    }
    
    if (best_solution_.feasible) {
        printf("\nBest cost: %.2f\n", best_solution_.total_cost);
        printf("Best hubs: ");
        for (size_t h = 0; h < best_solution_.hubs.size(); h++) {
            printf("%d ", best_solution_.hubs[h] + 1);
        }
        printf("\n");
    }
    
    solution = best_binary_solution_;
    return best_solution_.feasible;
}

bool RKOOptimizer::solveAllAlgorithms(const THLPData& data) {
    return solveTHLP(data, best_binary_solution_);
}

void RKOOptimizer::printSummary() const {
    printf("\n=== Summary ===\n");
    if (best_solution_.feasible) {
        printf("Best cost: %.2f\n", best_solution_.total_cost);
        printf("Hubs: ");
        for (size_t h = 0; h < best_solution_.hubs.size(); h++) {
            printf("%d ", best_solution_.hubs[h] + 1);
        }
        printf("\n");
    }
}

// ============================================================================
// 21. CAB PARSER
// ============================================================================

static bool parseCABFile(const std::string& filename, THLPData& data, 
                         int p = 3, double alpha = 0.8) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        printf("Error: Cannot open %s\n", filename.c_str());
        return false;
    }
    
    data.p = p;
    data.alpha = alpha;
    
    std::string line;
    std::vector<std::tuple<int, int, double, double> > entries;
    bool inData = false;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t first = line.find_first_not_of(" \t");
        if (first != std::string::npos && line[first] == '#') continue;
        if (line.find("param:") != std::string::npos) { inData = true; continue; }
        if (!inData) continue;
        if (line.find(";") != std::string::npos) continue;
        
        std::istringstream iss(line);
        int i, j;
        double W, C;
        if (iss >> i >> j >> W >> C) {
            entries.push_back(std::make_tuple(i-1, j-1, W, C));
        }
    }
    file.close();
    
    if (entries.empty()) {
        printf("Error: No data parsed\n");
        return false;
    }
    
    int maxNode = 0;
    for (std::vector<std::tuple<int, int, double, double> >::const_iterator it = entries.begin();
         it != entries.end(); ++it) {
        maxNode = std::max(maxNode, std::max(std::get<0>(*it), std::get<1>(*it)));
    }
    data.n = maxNode + 1;
    data.init(data.n);
    
    for (std::vector<std::tuple<int, int, double, double> >::const_iterator it = entries.begin();
         it != entries.end(); ++it) {
        int i = std::get<0>(*it);
        int j = std::get<1>(*it);
        data.W[i][j] = std::get<2>(*it);
        data.C[i][j] = std::get<3>(*it);
    }
    
    data.precompute();
    printf("Loaded CAB: %d nodes, p=%d, alpha=%.2f\n", data.n, p, alpha);
    return true;
}

// ============================================================================
// 22. ENTRY POINTS
// ============================================================================

bool solveTHLP(const THLPData& data, std::vector<double>& solution) {
    RKOConfig config;
    RKOOptimizer optimizer(config);
    return optimizer.solveTHLP(data, solution);
}

bool testAllAlgorithmsOnTHLP() {
    printf("\n============================================================\n");
    printf("   Testing All Algorithms on THLP\n");
    printf("============================================================\n");
    
    THLPData data;
    if (!parseCABFile("cab25.txt", data, 5, 0.5)) {
        return false;
    }
    
    RKOConfig config;
    config.num_algorithms = 8;
    config.num_runs = 4;
    config.enable_ipr = true;
    
    RKOOptimizer optimizer(config);
    return optimizer.solveAllAlgorithms(data);
}

// ============================================================================
// 23. MAIN
// ============================================================================

int main() {
    bool success = testAllAlgorithmsOnTHLP();
    printf("\n============================================================\n");
    printf("   Done. Success: %s\n", success ? "Yes" : "No");
    printf("============================================================\n");
    return success ? 0 : 1;
}