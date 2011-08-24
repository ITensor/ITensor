#define THIS_IS_MAIN
#include "hams.h"
#include <boost/python.hpp>
using namespace boost::python;

//Wrappers:
IndexVal Index_call(const Index& self, int i) { return self(i); }

void Matrix_set(Matrix& self, int i, int j, Real val) { self(i,j) = val; }

void ITensor_set1(ITensor& self, 
                  const IndexVal& iv1, 
                  Real val) 
                  { self(iv1) = val; }

void ITensor_set2(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  Real val) 
                  { self(iv1,iv2) = val; }

void ITensor_set3(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  Real val) 
                  { self(iv1,iv2,iv3) = val; }

void ITensor_set4(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  const IndexVal& iv4, 
                  Real val) 
                  { self(iv1,iv2,iv3,iv4) = val; }

void ITensor_set5(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  const IndexVal& iv4, 
                  const IndexVal& iv5, 
                  Real val) 
                  { self(iv1,iv2,iv3,iv4,iv5) = val; }

void ITensor_set6(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  const IndexVal& iv4, 
                  const IndexVal& iv5, 
                  const IndexVal& iv6, 
                  Real val) 
                  { self(iv1,iv2,iv3,iv4,iv5,iv6) = val; }

void ITensor_set7(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  const IndexVal& iv4, 
                  const IndexVal& iv5, 
                  const IndexVal& iv6, 
                  const IndexVal& iv7, 
                  Real val) 
                  { self(iv1,iv2,iv3,iv4,iv5,iv6,iv7) = val; }

void ITensor_set8(ITensor& self, 
                  const IndexVal& iv1, 
                  const IndexVal& iv2, 
                  const IndexVal& iv3, 
                  const IndexVal& iv4, 
                  const IndexVal& iv5, 
                  const IndexVal& iv6, 
                  const IndexVal& iv7, 
                  const IndexVal& iv8, 
                  Real val) 
                  { self(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8) = val; }

void InitState_set(InitState& self, int j, const IQIndexVal& iv) { self(j) = iv; }

class SpinOneModel
{
    SpinOne::Model model_;
public:
    const SpinOne::Model& model() const { return model_; }

    SpinOneModel(int N) : model_(N) { }

    int NN() const { return model_.NN(); }

    IQIndexVal Up(int j) const { return model_.Up(j); }
    IQIndexVal Z0(int j) const { return model_.Z0(j); }
    IQIndexVal Dn(int j) const { return model_.Dn(j); }

    ITensor id(int j) const { return model_.id(j); }
};

class MPSWrapper
{
    const BaseModel& model;
    MPS mps;
public:
    MPSWrapper(const SpinOneModel& s1m) : model(s1m.model()),mps(model) { }
    MPSWrapper(const SpinOneModel& s1m, const InitState& initState) : model(s1m.model()),mps(model,initState) { }

    int NN() const { return model.NN(); }
    ITensor A(int j) const { return mps.AA(j); }

    void position(int i) { mps.position(i); }

    ITensor projectOp(int j, Direction dir, const ITensor& P, const ITensor& Op)
    { ITensor res; mps.projectOp(j,dir,P,Op,res); return res; }

    Real bondDavidson(int b, const ITensor& mpoh, const ITensor& LH, const ITensor& RH, int niter, int debuglevel, Direction dir)
    { mps.bondDavidson(b,mpoh,LH,RH,niter,debuglevel,dir); }

    friend inline ostream& operator<<(ostream& s, const MPSWrapper& M) { s << M.mps; return s; }
};

class MPOWrapper
{
    const BaseModel& model;
    MPO mpo;
public:
    MPOWrapper(const MPO& mpo_) : model(mpo_.model()),mpo(mpo_) { }
    MPOWrapper(const SpinOneModel& s1m) : model(s1m.model()),mpo(model) { }

    int NN() const { return model.NN(); }
    ITensor A(int j) const { return mpo.AA(j); }

    friend inline ostream& operator<<(ostream& s, const MPOWrapper& M) { s << M.mpo; return s; }
};

MPOWrapper SpinOneHeisenberg(const SpinOneModel& s1m) { MPO H = SpinOne::Heisenberg(s1m.model())(); MPOWrapper res(H); return res; }


BOOST_PYTHON_MODULE(itensor)
{
    class_<Matrix>("MatrixBase", init<int,int>())
    .def("set",&Matrix_set)
    .def(self_ns::str(self))
    ;

    class_<Index>("Index", init<string,int>())
    .add_property("m", &Index::m)
    .add_property("name", &Index::name)
    .add_property("is_null", &Index::is_null)
    .def("__call__",&Index_call)
    .def(self_ns::str(self))
    ;

    class_<IndexVal>("IndexVal", init<Index,int>())
    .def_readonly("ind",&IndexVal::ind)
    .def_readonly("i",&IndexVal::i)
    .def(self_ns::str(self))
    ;

    class_<ITensor>("ITensor")
    .def(init<Index>())
    .def(init<Index,Index>())
    .def(init<Index,Index,Matrix>())
    .def(init<Index,Index,Index>())
    .def(init<Index,Index,Index,Index>())
    .def("set",&ITensor_set1).def("set",&ITensor_set2)
    .def("set",&ITensor_set3).def("set",&ITensor_set4)
    .def("set",&ITensor_set5).def("set",&ITensor_set6)
    .def("set",&ITensor_set7).def("set",&ITensor_set8)
    .def("__len__",&ITensor::vec_size)
    .def("norm",&ITensor::norm)
    .def(self_ns::str(self))
    .def(self * self)
    .def(self *= self)
    ;

    enum_<Direction>("Direction")
    .value("Fromleft",Fromleft)
    .value("Fromright",Fromright)
    ;

    class_<IQIndexVal>("IQIndexVal", init<IQIndex,int>())
    .def(self_ns::str(self))
    ;

    class_<SpinOneModel>("SpinOneModel",init<int>())
    .def("Up",&SpinOneModel::Up)
    .def("Z0",&SpinOneModel::Z0)
    .def("Dn",&SpinOneModel::Dn)
    ;

    class_<InitState>("InitState",init<int>())
    .def("set",&InitState_set)
    ;

    class_<MPSWrapper>("MPS",init<SpinOneModel>())
    .def(init<SpinOneModel,InitState>())
    .def("A",&MPSWrapper::A)
    .def("position",&MPSWrapper::position)
    .def("projectOp",&MPSWrapper::projectOp)
    .def("bondDavidson",&MPSWrapper::bondDavidson)
    .def("__len__",&MPSWrapper::NN)
    .def(self_ns::str(self))
    ;

    class_<MPOWrapper>("MPO",init<SpinOneModel>())
    .def("A",&MPOWrapper::A)
    .def("__len__",&MPOWrapper::NN)
    .def(self_ns::str(self))
    ;

    def("SpinOneHeisenberg",&SpinOneHeisenberg);
}
