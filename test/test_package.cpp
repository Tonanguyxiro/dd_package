#include <memory>

#include "DDpackage.h"
#include "DDexport.h"
#include "util.h"
#include "gtest/gtest.h"
#include <regex>

TEST(DDPackageTest, TrivialTest) {
    auto dd = std::make_unique<dd::Package>();

    short line[2] = {2};
    dd::Edge x_gate = dd->makeGateDD(Xmat, 1, line);
    dd::Edge h_gate = dd->makeGateDD(Hmat, 1, line);

    ASSERT_EQ(dd->getValueByPath(h_gate, "0"), (dd::ComplexValue{dd::SQRT_2, 0}));

    dd::Edge zero_state = dd->makeZeroState(1);
    dd::Edge h_state = dd->multiply(h_gate, zero_state);
    dd::Edge one_state = dd->multiply(x_gate, zero_state);

    ASSERT_EQ(dd->fidelity(zero_state, one_state), 0.0);
    ASSERT_NEAR(dd->fidelity(zero_state, h_state), 0.5, dd::ComplexNumbers::TOLERANCE);
    ASSERT_NEAR(dd->fidelity(one_state, h_state), 0.5, dd::ComplexNumbers::TOLERANCE);
}

TEST(DDPackageTest, BellState) {
    auto dd = std::make_unique<dd::Package>();

    short line[2] = {-1,2};
    dd::Edge h_gate = dd->makeGateDD(Hmat, 2, line);
    dd::Edge cx_gate = dd->makeGateDD({Xmat[0][0], Xmat[0][1], Xmat[1][0], Xmat[1][1]}, 2, {2,1});
    dd::Edge zero_state = dd->makeZeroState(2);

    dd::Edge bell_state = dd->multiply(dd->multiply(cx_gate, h_gate), zero_state);

    ASSERT_EQ(dd->getValueByPath(bell_state, "00"), (dd::ComplexValue{dd::SQRT_2, 0}));
    ASSERT_EQ(dd->getValueByPath(bell_state, "02"), (dd::ComplexValue{0, 0}));
    ASSERT_EQ(dd->getValueByPath(bell_state, "20"), (dd::ComplexValue{0, 0}));
    ASSERT_EQ(dd->getValueByPath(bell_state, "22"), (dd::ComplexValue{dd::SQRT_2, 0}));

    ASSERT_DOUBLE_EQ(dd->fidelity(zero_state, bell_state), 0.5);
    dd->printDD(bell_state, 64);
}

TEST(DDPackageTest, IdentifyTrace) {
    auto dd = std::make_unique<dd::Package>();
    auto fullTrace = dd->trace(dd->makeIdent(0, 3));

    ASSERT_EQ(fullTrace, (dd::ComplexValue{16,0}));
}

TEST(DDPackageTest, BellStateAmplitudeComparison) {
    auto dd = std::make_unique<dd::Package>();

    short line[2] = {-1,2};
    dd::Edge h_gate = dd->makeGateDD(Hmat, 2, line);
    dd::Edge cx_gate = dd->makeGateDD({Xmat[0][0], Xmat[0][1], Xmat[1][0], Xmat[1][1]}, 2, {2,1});
    dd::Edge zero_state = dd->makeZeroState(2);

    dd::Edge bell_state = dd->multiply(dd->multiply(cx_gate, h_gate), zero_state);

    dd->setMode(dd::Mode::Matrix);

    dd::Edge bell_state_matrix = dd->multiply(dd->multiply(cx_gate, h_gate), zero_state);
        
    EXPECT_FALSE(dd->equals(bell_state, bell_state_matrix));

    unsigned short nqubits = 2;
    std::map<std::string, dd::ComplexValue> result_amplitudes;
    std::map<std::string, dd::ComplexValue> ref_amplitudes;

    std::string elements(nqubits, '0');
    dd->getAllAmplitudes(bell_state, result_amplitudes, nqubits - 1, elements);
    dd->getAllAmplitudes(bell_state_matrix, ref_amplitudes, nqubits - 1, elements);
    EXPECT_TRUE(dd->compareAmplitudes(ref_amplitudes, result_amplitudes, false));               
}

TEST(DDPackageTest, BellStateSerialization) {
    auto dd = std::make_unique<dd::Package>();

    short line[2] = {-1,2};
    dd::Edge h_gate = dd->makeGateDD(Hmat, 2, line);
    dd::Edge cx_gate = dd->makeGateDD({Xmat[0][0], Xmat[0][1], Xmat[1][0], Xmat[1][1]}, 2, {2,1});
    dd::Edge zero_state = dd->makeZeroState(2);

    dd::Edge bell_state = dd->multiply(dd->multiply(cx_gate, h_gate), zero_state);

    dd::serialize(bell_state, "bell_state", true);
    // dd->export2Dot(bell_state, "bell_state_graph", true);

    dd::Edge result = dd::deserialize(dd, "bell_state");

    // dd->export2Dot(result, "bell_state_graph_deserialized", true);

    EXPECT_TRUE(dd->equals(bell_state, result));
}

TEST(DDPackageTest, DeleteFirstEdge) {
    auto dd = std::make_unique<dd::Package>();

    dd::Edge zero = dd->makeZeroState(3);
    /*
    dd::Edge zero_deleted = dd->deleteEdge(zero, 2, 0);
    EXPECT_TRUE(dd->equals(zero_deleted, dd::Package::DDzero));        
    EXPECT_FALSE(dd->equals(zero, zero_deleted));        
    */    
    short line[3] = {-1, -1, 2};
    dd::Edge h_gate = dd->makeGateDD(Hmat, 3, line);
    dd::Edge hzero = dd->multiply(h_gate, zero);    
    dd::Edge hzero_deleted = dd->deleteEdge(hzero, 2, 0);
    
    EXPECT_TRUE(dd->equals(hzero_deleted.p->e[0], dd::Package::DDzero));    
}