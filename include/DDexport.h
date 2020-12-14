/*
 * This file is part of IIC-JKU DD package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDexport_H
#define DDexport_H

#include "DDpackage.h"

#include <iomanip>

namespace dd {
	constexpr char SERIALIZATION_VERSION[] = "0.1";

	struct RGB {
		fp R=0., G=0., B=0.;
		RGB(fp R, fp G, fp B): R(R), G(G), B(B) {};

		std::ostream& printHex(std::ostream& os) const {
			std::ostringstream oss{};
			oss.flags(std::ios_base::hex);
			oss.fill('0');
			oss << std::setw(2) << short(R*255) << std::setw(2) << short(G*255) << std::setw(2) << short(B*255);
			os << oss.str();
			return os;
		}

		friend std::ostream& operator<<(std::ostream& os, const RGB& rgb) {
			return rgb.printHex(os);
		}
	};

	fp hueToRGB(fp hue);
	RGB hlsToRGB(const fp& h, const fp& l, const fp& s);

	RGB colorFromPhase(const Complex& a);
	fp thicknessFromMagnitude (const Complex& a);

	std::ostream& header(const Edge& e, std::ostream& os, bool edgeLabels, bool colored = false);

	std::ostream& matrixNodeMatrixAndXlabel(const Edge& e, std::ostream& os);
	std::ostream& matrixNodeMiddleVar(const Edge& e, std::ostream& os);
	std::ostream& classicMatrixNode(const Edge& e, std::ostream& os);

	std::ostream& vectorNode(const Edge& e, std::ostream& os);
	std::ostream& vectorNodeVectorLook(const Edge& e, std::ostream& os);
	std::ostream& classicVectorNode(const Edge& e, std::ostream& os);

	std::ostream& matrixEdge(const Edge& from, const Edge& to, short idx, std::ostream& os, bool edgeLabels=false, bool classic=false, bool colored=false);
	std::ostream& coloredMatrixEdge(const Edge& from, const Edge& to, short idx, std::ostream& os, bool edgeLabels=false, bool classic=false);

	std::ostream& vectorEdge(const Edge& from, const Edge& to, short idx, std::ostream& os, bool edgeLabels=false, bool classic=false, bool colored=false);

	void toDot(const Edge& e, std::ostream& os, bool isVector = false, bool colored=true, bool edgeLabels=false, bool classic=false);
	void export2Dot(Edge basic, const std::string& outputFilename, bool isVector = false, bool colored=true, bool edgeLabels=false, bool classic=false, bool show = true);

	ComplexValue toComplexValue(const std::string& real_str, std::string imag_str);

	void serialize(Edge basic, const std::string& outputFilename, bool isVector = false);
	void serialize(Edge basic, std::ostream& oss, bool isVector = false);
	dd::Edge deserialize(std::unique_ptr<dd::Package>& dd, const std::string& inputFilename);
	dd::Edge deserialize(std::unique_ptr<dd::Package>& dd, std::istream& ifs);

	void exportAmplitudes(std::unique_ptr<dd::Package>& dd, Edge basic, const std::string& outputFilename);
	void exportAmplitudesRec(const Edge& node, std::ostream& oss, std::string path, Complex& amplitude);
	void exportAmplitudes(std::unique_ptr<dd::Package>& dd, Edge basic, std::ostream& oss);
}


#endif //DDexport_H
