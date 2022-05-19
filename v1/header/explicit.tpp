#ifndef QSL_EXPLICIT_TPP
#define QSL_EXPLICIT_TPP

#include "qsl.hpp"

// List of all explicit instantiations

template class qsl::basic<float, true, qsl::seq>;
template class qsl::basic<float, true, qsl::omp>;
template class qsl::basic<float, true, qsl::opt>;
template class qsl::basic<float, false, qsl::seq>;
template class qsl::basic<float, false, qsl::omp>;
template class qsl::basic<float, false, qsl::opt>;
template class qsl::basic<double, true, qsl::seq>;
template class qsl::basic<double, true, qsl::omp>;
template class qsl::basic<double, true, qsl::opt>;
template class qsl::basic<double, false, qsl::seq>;
template class qsl::basic<double, false, qsl::omp>;
template class qsl::basic<double, false, qsl::opt>;
template class qsl::basic<long double, true, qsl::seq>;
template class qsl::basic<long double, true, qsl::omp>;
template class qsl::basic<long double, true, qsl::opt>;
template class qsl::basic<long double, false, qsl::seq>;
template class qsl::basic<long double, false, qsl::omp>;
template class qsl::basic<long double, false, qsl::opt>;
template class qsl::resize<float, true, qsl::seq>;
template class qsl::resize<float, true, qsl::omp>;
template class qsl::resize<float, true, qsl::opt>;
template class qsl::resize<float, false, qsl::seq>;
template class qsl::resize<float, false, qsl::omp>;
template class qsl::resize<float, false, qsl::opt>;
template class qsl::resize<double, true, qsl::seq>;
template class qsl::resize<double, true, qsl::omp>;
template class qsl::resize<double, true, qsl::opt>;
template class qsl::resize<double, false, qsl::seq>;
template class qsl::resize<double, false, qsl::omp>;
template class qsl::resize<double, false, qsl::opt>;
template class qsl::resize<long double, true, qsl::seq>;
template class qsl::resize<long double, true, qsl::omp>;
template class qsl::resize<long double, true, qsl::opt>;
template class qsl::resize<long double, false, qsl::seq>;
template class qsl::resize<long double, false, qsl::omp>;
template class qsl::resize<long double, false, qsl::opt>;
template class qsl::number<float, true, qsl::seq>;
template class qsl::number<float, true, qsl::omp>;
template class qsl::number<float, true, qsl::opt>;
template class qsl::number<float, false, qsl::seq>;
template class qsl::number<float, false, qsl::omp>;
template class qsl::number<float, false, qsl::opt>;
template class qsl::number<double, true, qsl::seq>;
template class qsl::number<double, true, qsl::omp>;
template class qsl::number<double, true, qsl::opt>;
template class qsl::number<double, false, qsl::seq>;
template class qsl::number<double, false, qsl::omp>;
template class qsl::number<double, false, qsl::opt>;
template class qsl::number<long double, true, qsl::seq>;
template class qsl::number<long double, true, qsl::omp>;
template class qsl::number<long double, true, qsl::opt>;
template class qsl::number<long double, false, qsl::seq>;
template class qsl::number<long double, false, qsl::omp>;
template class qsl::number<long double, false, qsl::opt>;

#endif
