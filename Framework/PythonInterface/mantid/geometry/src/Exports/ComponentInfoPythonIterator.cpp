// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidPythonInterface/api/ComponentInfoPythonIterator.h"

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/module.hpp>

using Mantid::PythonInterface::ComponentInfoPythonIterator;
using namespace boost::python;

// Export ComponentInfoPythonIterator
void export_ComponentInfoPythonIterator() {

  // Export to Python
  class_<ComponentInfoPythonIterator>("ComponentInfoPythonIterator", no_init)
      .def("__iter__", objects::identity_function())
#ifdef IS_PY3K
      .def("__next__", &ComponentInfoPythonIterator::next)
#else
      .def("next", &ComponentInfoPythonIterator::next,
           return_value_policy<copy_const_reference>())
#endif
      ;
  /*
   Return value policy for next is to copy the const reference. Copy by value is
   essential for python 2.0 compatibility because items (DetectorInfoItem) will
   outlive their iterators if declared as part of for loops. There is no good
   way to deal with this other than to force a copy so that internals of the
   item are not also corrupted. Future developers may wish to choose a separte
   policy for python 3.0 where this is not a concern, and const ref returns
   would be faster.
  */
}
