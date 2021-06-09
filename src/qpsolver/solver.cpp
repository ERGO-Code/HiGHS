#include "solver.hpp"

#include <algorithm>
#include <map>

#include "lp_data/HighsAnalysis.h"
#include "Highs.h"

#include "basis.hpp"
#include "feasibility/feasibility_highs.hpp"
#include "gradient.hpp"
#include "redhes/factor.hpp"
#include "pricing/dantzigpricing.hpp"
#include "pricing/devexpricing.hpp"
#include "pricing/devexharrispricing.hpp"
#include "pricing/steepestedgepricing.hpp"
#include "ratiotest/ratiotest.hpp"
#include "reducedgradient.hpp"
#include "snippets.hpp"
#include "reducedcosts.hpp"

#include "instance.hpp"

#include "nullspace.hpp"


void Solver::solve() {
	CrashSolution* crash;
	computestartingpoint(runtime, crash);   
   Basis basis(runtime, crash->active, crash->rowstatus, crash->inactive);
   solve(crash->primal, crash->rowact,  basis);
}

Solver::Solver(Runtime& rt) : runtime(rt) {

}

void Solver::loginformation(Runtime& rt, Basis& basis, Nullspace& ns, NewCholeskyFactor& factor) {
	rt.statistics.iteration.push_back(rt.statistics.num_iterations);
	rt.statistics.nullspacedimension.push_back(rt.instance.num_var - basis.getnumactive());
	rt.statistics.objval.push_back(rt.instance.objval(rt.primal));
	rt.statistics.time.push_back(std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - rt.statistics.time_start).count());
	SumNum sm = rt.instance.sumnumprimalinfeasibilities(rt.primal, rt.rowactivity);
	rt.statistics.sum_primal_infeasibilities.push_back(sm.sum);	
	rt.statistics.num_primal_infeasibilities.push_back(sm.num);
	rt.statistics.density_factor.push_back(factor.density());
	rt.statistics.density_nullspace.push_back(ns.density());
}

void tidyup(Vector& p, Vector& rowmove, Basis& basis, Runtime& runtime) {
	for (unsigned acon : basis.getactive()) {
		if (acon >= runtime.instance.num_con) {
			p.value[acon - runtime.instance.num_con] = 0.0;
		} else {
			rowmove.value[acon] = 0.0;
		}
	}
}

void recomputexatfsep(Runtime& runtime) {

}

void computerowmove(Runtime& runtime, Basis& basis, Vector& p, Vector& rowmove) {
	runtime.instance.A.mat_vec(p, rowmove);
	return;
	// rowmove.reset();
	MatrixBase& Atran = runtime.instance.A.t();
	Atran.vec_mat(p, rowmove);
	return;
	for (int i=0; i<runtime.instance.num_con; i++) {
		if (basis.getstatus(i) == BasisStatus::Default) {
			// check with assertions, is it really the same?
			double val = p.dot(&Atran.index[Atran.start[i]], &Atran.value[Atran.start[i]], Atran.start[i+1] - Atran.start[i]);
			// Vector col = Atran.extractcol(i);
			// val = col * p;
			
			// assert(rowmove.value[i] == val);
			rowmove.value[i] = val;
		} else {
			rowmove.value[i] = 0;
		}
	}
	rowmove.resparsify();
} 

// VECTOR
Vector& computesearchdirection_minor(Runtime& rt, Nullspace& ns, Basis& bas,  NewCholeskyFactor& cf, ReducedGradient& redgrad, Vector& p) {		
	Vector g2 = -redgrad.get(); 
	g2.sanitize();
	cf.solve(g2);

	g2.sanitize();

	return ns.Zprod(g2, p);
}

// VECTOR
Vector& computesearchdirection_major(Runtime& runtime, Nullspace& ns, Basis& basis, NewCholeskyFactor& factor, Vector& yp, Gradient& gradient, Vector& gyp, Vector& l, Vector& p) {		
	runtime.instance.Q.mat_vec(yp, gyp);
	if (basis.getnumactive() < runtime.instance.num_var) {
		ns.Ztprod(gyp, l);
		factor.solveL(l);
		Vector v = l;
		factor.solveLT(v);
		ns.Zprod(v, p);
		return p.saxpy(-1.0, 1.0, yp);
	} else {
		return p.repopulate(yp).scale(-gradient.getGradient().dot(yp));
		// return -yp; 
	}
}

void Solver::solve(const Vector& x0, const Vector& ra, Basis& b0) {
	runtime.statistics.time_start = std::chrono::high_resolution_clock::now();
	Basis& basis = b0;
	runtime.primal = x0;

	// TODO: remove redundant equations before starting 
	// HOWTO: from crash start, check all (near-)equality constraints (not bounds). if the residual is 0 (or near-zero?), remove constraint
	Nullspace ns(runtime, basis);
	Gradient gradient(runtime);
	ReducedCosts redcosts(runtime, basis, gradient);
	ReducedGradient redgrad(runtime, ns, gradient);
	NewCholeskyFactor factor(runtime, basis, ns);
	runtime.instance.A.mat_vec(runtime.primal, runtime.rowactivity);
	std::unique_ptr<Pricing> pricing = std::unique_ptr<Pricing>(new DevexPricing(runtime, basis, redcosts));
	
	Vector p(runtime.instance.num_var);
	Vector rowmove(runtime.instance.num_con);

	Vector buffer_yp(runtime.instance.num_var);
	Vector buffer_gyp(runtime.instance.num_var);
	Vector buffer_l(runtime.instance.num_var);

	Vector buffer_Qp(runtime.instance.num_var);

	bool atfsep = basis.getnumactive() == runtime.instance.num_var;
	while (true && runtime.statistics.num_iterations < runtime.settings.iterationlimit
	&& std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - runtime.statistics.time_start).count() < runtime.settings.timelimit) {
		// LOGGING
		if (runtime.statistics.num_iterations % runtime.settings.reportingfequency == 0) {
			loginformation(runtime, basis, ns, factor);
			runtime.endofiterationevent.fire(runtime);
		}
		runtime.statistics.num_iterations++;

		double maxsteplength = 1.0;
		if (atfsep) {
			int minidx = pricing->price(runtime.primal, gradient.getGradient());
			// printf("%u -> ", minidx);
			if (minidx == -1) {
				runtime.status = ProblemStatus::OPTIMAL;
				break;
			}
			
			ns.expand_computenewcol(minidx, buffer_yp);
			buffer_l.dim = basis.getnuminactive();
			computesearchdirection_major(runtime, ns, basis, factor, buffer_yp, gradient, buffer_gyp, buffer_l, p);
			basis.deactivate(minidx);
			computerowmove(runtime, basis, p, rowmove);
			tidyup(p, rowmove, basis, runtime);
			maxsteplength = std::numeric_limits<double>::infinity();
			if (runtime.instance.Q.mat.value.size() > 0) {
				double denominator = p * runtime.instance.Q.mat_vec(p, buffer_Qp);
				if (fabs(denominator) > 10E-5) {
					double numerator = -(p * gradient.getGradient());
					if (numerator < 0.0) {
						// printf("numerator < 0: %lf\n", numerator);
						maxsteplength = 0.0;
					} else {
						maxsteplength = numerator / denominator;
					}
				}
				factor.expand(buffer_yp, buffer_gyp, buffer_l);
			}
			redgrad.expand(buffer_yp);
			ns.expand_appendnewcol(buffer_yp);
		} else {
			computesearchdirection_minor(runtime, ns, basis, factor, redgrad, p);
			computerowmove(runtime, basis, p, rowmove);
			tidyup(p, rowmove, basis, runtime);
		}
		
				
		if (p.norm2() < runtime.settings.pnorm_zero_threshold || maxsteplength == 0.0) {
			atfsep = true;
		} else {
			RatiotestResult stepres = runtime.settings.ratiotest->ratiotest(runtime.primal, p, runtime.rowactivity, rowmove, runtime.instance, maxsteplength);
			// printf("%u, alpha= %lf\n", stepres.limitingconstraint,stepres.alpha);
			if (stepres.limitingconstraint != -1) {
				// Vector d = computed(runtime, ns, basis, stepres.limitingconstraint);
				NullspaceReductionResult nrr = ns.reduce(runtime, stepres.limitingconstraint);
				if (runtime.instance.Q.mat.value.size() > 0) {
					factor.reduce(nrr);
				}
				redgrad.reduce(nrr);
				redgrad.update(stepres.alpha, false);	
				
				
				basis.activate(runtime, stepres.limitingconstraint, stepres.nowactiveatlower ? BasisStatus::ActiveAtLower : BasisStatus::ActiveAtUpper, nrr.constrainttodrop, pricing.get());
				if (basis.getnumactive() != runtime.instance.num_var) {
					atfsep = false;
				}
			} else {
				if (stepres.limitingconstraint == std::numeric_limits<double>::infinity()) {
					// unbounded
					runtime.status = ProblemStatus::UNBOUNDED;
				}
				atfsep = true;
				redgrad.update(stepres.alpha, false);
			}
			
			gradient.update(p, stepres.alpha);
			redcosts.update();
			
			runtime.primal.saxpy(stepres.alpha, p);
			runtime.rowactivity.saxpy(stepres.alpha, rowmove);
		}

	}

	loginformation(runtime, basis, ns, factor);
	runtime.endofiterationevent.fire(runtime);

	Vector lambda =redcosts.getReducedCosts();
	for (auto e: basis.getactive()) {
		int indexinbasis = basis.getindexinfactor()[e];
		if (e >= runtime.instance.num_con) {
			// active variable bound
			int var = e - runtime.instance.num_con;
			runtime.dualvar.value[var] = lambda.value[indexinbasis];
		} else {
			runtime.dualcon.value[e] = lambda.value[indexinbasis];
		}
	}

	if (basis.getnumactive() == runtime.instance.num_var) {
		runtime.primal = basis.recomputex(runtime.instance);
	}
	// x.report("x");
	runtime.statistics.time_end = std::chrono::high_resolution_clock::now();
}