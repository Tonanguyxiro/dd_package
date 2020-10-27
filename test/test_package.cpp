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

TEST(DDPackageTest, IdentityTrace) {
    auto dd = std::make_unique<dd::Package>();
    auto fullTrace = dd->trace(dd->makeIdent(0, 3));

    ASSERT_EQ(fullTrace, (dd::ComplexValue{16,0}));
}

TEST(DDPackageTest, StateGenerationManipulation) {
	auto dd = std::make_unique<dd::Package>();

	auto b = std::bitset<dd::MAXN>{2};
	auto e = dd->makeBasisState(6, b);
	auto f = dd->makeBasisState(6, {dd::BasisStates::zero,
								            dd::BasisStates::one,
								            dd::BasisStates::plus,
								            dd::BasisStates::minus,
								            dd::BasisStates::left,
								            dd::BasisStates::right});
	dd->incRef(e);
	dd->incRef(f);
	dd->incRef(e);
	auto g = dd->add(e, f);
	auto h = dd->transpose(g);
	auto i = dd->conjugateTranspose(f);
	dd->decRef(e);
	dd->decRef(f);
	auto j = dd->kronecker(h, i);
	dd->incRef(j);
	dd->printActive(6);
	dd->printUniqueTable(6);
	dd->printInformation();
}

TEST(DDPackageTest, BellStateSerialization) {
    auto dd = std::make_unique<dd::Package>();

    short line[2] = {-1,2};
    dd::Edge h_gate = dd->makeGateDD(Hmat, 2, line);
    dd::Edge cx_gate = dd->makeGateDD({Xmat[0][0], Xmat[0][1], Xmat[1][0], Xmat[1][1]}, 2, {2,1});
    dd::Edge zero_state = dd->makeZeroState(2);

    dd::Edge bell_state = dd->multiply(dd->multiply(cx_gate, h_gate), zero_state);

    dd::serialize(bell_state, "bell_state", true);
    dd->export2Dot(bell_state, "bell_state_graph", true);

    dd::Edge result = dd::deserialize(dd, "bell_state");

    dd->export2Dot(result, "bell_state_graph_deserialized", true);

    EXPECT_TRUE(dd->equals(bell_state, result));

    /*
    std::string complex_real_regex = "([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)?";
    std::string complex_imag_regex = "([+-]?(?:(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)?[iI])?";

    std::string test = "0 0 (-1 1) () () ()";
    // std::string test = "2 1 (0 0.8+5i) () () ()";

    std::string edge_regex = " \\(((-?\\d+) (" + complex_real_regex + " ?" + complex_imag_regex + "))?\\)";
    // std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + "){4} *#?.*"); // TODO {4} overwrites groups
    std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ") *#?.*");
    std::smatch m;

    if(std::regex_match(test, m, e, std::regex_constants::match_prev_avail)) {

        // match 1: node_idx 
        // match 2: qubit_idx 

        // repeats for every edge
        // match 3: edge content
        // match 4: edge_target_idx 
        // match 5: real + imag
        // match 6: real
        // match 7: imag




        std::cout << "MATCH" << std::endl;
        for(int i = 0; i < m.size(); i++) {
        std::cout << i << ": " << m.str(i) << std::endl;

        }
    }
    */
   
   // http://www.cplusplus.com/reference/regex/ECMAScript/
    // std::string test = "2 1 (0 0.7071067811865476) () (1 0.7071067811865476) ()";
    /*
    std::string test = "2 1 (0 0.8)";
    std::string edge_regex = " \\(\\d+ 0.\\d+\\)";
    std::regex e ("(\\d+) (\\d+)" + edge_regex);
    */
   /*
    std::string test = "2 1 (0 0.8) () () () # ";
    std::string edge_regex = " \\(((\\d+) (0.\\d+))?\\)";
    std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + "){4} *#.*");
    std::smatch m;

    if(std::regex_match(test, m, e)) {
        std::cout << "MATCH" << std::endl;
        std::cout << m.str(1) << std::endl;
        std::cout << m.str(2) << std::endl;
        std::cout << m.str(3) << std::endl;
        std::cout << m.str(4) << std::endl;
        std::cout << m.str(5) << std::endl;
        std::cout << m.str(6) << std::endl;
        std::cout << m.str(7) << std::endl;
    }
    */
    
    /*
    std::string complex_real_regex = "([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)?";
    std::string complex_imag_regex = "([+-]?(?:(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)?[iI])?";

    std::string test = "2 1 (0 0.8+5i) (1 5+5i) () () # ";
    std::string edge_regex = " \\(((\\d+) (" + complex_real_regex + " ?" + complex_imag_regex + "))?\\)";
    // std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + "){4} *#.*"); // TODO {4} overwrites groups
    std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ") *#.*");
    std::smatch m;

    if(std::regex_match(test, m, e, std::regex_constants::match_prev_avail)) {

        // match 1: node_idx 
        // match 2: qubit_idx 

        // repeats for every edge
        // match 3: edge content
        // match 4: edge_target_idx 
        // match 5: real + imag
        // match 6: real
        // match 7: imag




        std::cout << "MATCH" << std::endl;
        for(int i = 0; i < m.size(); i++) {
        std::cout << i << ": " << m.str(i) << std::endl;

        }
    }
    */
}
