#include "parallel/HighsParallel.h"
#include "lp_data/HighsCallback.h"
#include "ipm/ipx/control.h"
#include <iostream>

namespace ipx {

Control::Control() {
    // When failbit is set, the stream evaluates to false.
    dummy_.setstate(std::ios::failbit);
}

Int Control::InterruptCheck(const Int ipm_iteration_count) const {
    HighsTaskExecutor::getThisWorkerDeque()->checkInterrupt();
    if (parameters_.time_limit >= 0.0 &&
        parameters_.time_limit < timer_.Elapsed())
        return IPX_ERROR_time_interrupt;
    // The pointer callback_ should not be null, since that indicates
    // that it's not been set
    assert(callback_);
    if (callback_) {
      if (callback_->user_callback && callback_->active[kCallbackIpmInterrupt]) {
	callback_->clearHighsCallbackDataOut();
	callback_->data_out.ipm_iteration_count = ipm_iteration_count;
	if (callback_->callbackAction(kCallbackIpmInterrupt,
				      "IPM interrupt"))
	  return IPX_ERROR_user_interrupt;
      }
    }
    return 0;
}

std::stringstream Control::hLoggingStream() const {
  std::stringstream logging;
  logging.str(std::string());
  return logging;
}

void Control::hLog(std::stringstream& logging) const {
  if (parameters_.highs_logging) {
    highsLogUser(parameters_.log_options, HighsLogType::kInfo, "%s", logging.str().c_str());
    output_ << "output_ << " << logging.str();
  } else {
    output_ << logging.str();
  }
  logging.str(std::string());
}

std::ostream& Control::Log() const {
  if (parameters_.highs_logging) {
    //     std::string myString = output_.str();
     return dummy_;
  } else {
    return output_;
  }
}

std::ostream& Control::IntervalLog() const {
    if (parameters_.print_interval >= 0.0 &&
        interval_.Elapsed() >= parameters_.print_interval) {
        interval_.Reset();
        return output_;
    } else {
        return dummy_;
    }
}

std::ostream& Control::Debug(Int level) const {
    if (parameters_.debug >= level)
        return output_;
    else
        return dummy_;
}

void Control::ResetPrintInterval() const {
    interval_.Reset();
}

double Control::Elapsed() const {
    return timer_.Elapsed();
}

const Parameters& Control::parameters() const {
    return parameters_;
}

void Control::parameters(const Parameters& new_parameters) {
    parameters_ = new_parameters;
    MakeStream();
}

void Control::callback(HighsCallback* callback) {
    callback_ = callback;
}

void Control::OpenLogfile() {
    logfile_.close();
    const char* filename = parameters_.logfile;
    if (filename && filename[0])
        logfile_.open(filename, std::ios_base::out | std::ios_base::app);
    MakeStream();
}

void Control::CloseLogfile() {
    logfile_.close();
    MakeStream();
}

void Control::ResetTimer() {
    timer_.Reset();
}

void Control::MakeStream() {
    output_.clear();
    if (parameters_.display)
        output_.add(std::cout);
    if (logfile_.is_open())
        output_.add(logfile_);
}

std::string Format(Int i, int width) {
    std::ostringstream s;
    s.width(width);
    s << i;
    return s.str();
}

std::string Format(const char* c, int width) {
    std::ostringstream s;
    s.width(width);
    s << c;
    return s.str();
}

std::string Format(double d, int width, int prec,
                   std::ios_base::fmtflags floatfield) {
    std::ostringstream s;
    s.precision(prec);
    s.width(width);
    s.setf(floatfield, std::ios_base::floatfield);
    s << d;
    return s.str();
}

}  // namespace ipx
