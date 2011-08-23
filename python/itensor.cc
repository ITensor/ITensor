#define THIS_IS_MAIN
#include "tensor.h"
#include <boost/python.hpp>
using namespace boost::python;

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ITensor_set_from_indexval_overloads, ITensor::set_from_indexval, 1, 8)

BOOST_PYTHON_MODULE(itensor)
{
    class_<Index>("Index", init<string,int>())
    .add_property("m", &Index::m)
    .add_property("name", &Index::name)
    .add_property("is_null", &Index::is_null)
    .def(self_ns::str(self))
    ;

    class_<IndexVal>("IndexVal", init<Index,int>())
    .def_readonly("ind",&IndexVal::ind)
    .def_readonly("i",&IndexVal::i)
    .def(self_ns::str(self))
    ;

    class_<ITensor>("ITensor", init<Index>())
    .def("__len__",&ITensor::vec_size)
    .def("norm",&ITensor::norm)
    //.def("set",&ITensor::set_from_indexval,ITensor_set_from_indexval_overloads(args("iv1","iv2","iv3","iv4","iv5","iv6","iv7","iv8"),"Docstring"),return_value_policy<copy_non_const_reference>())
    //.def("set",&ITensor::from_one_iv)
    //.def("set_ref",&ITensor::from_one_iv_ref,return_value_policy<copy_non_const_reference>())
    .def(self_ns::str(self))
    .def(self * self)
    ;
}
