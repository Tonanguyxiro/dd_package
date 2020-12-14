/*
 * This file is part of IIC-JKU DD package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#include "DDexport.h"

#include <stack>
#include <regex>

namespace dd {
	std::ostream& header(const Edge& e, std::ostream& os, bool edgeLabels, bool colored) {
		os << "digraph \"DD\" {graph[];node[shape=plain];edge[arrowhead=none]\n";
		os << "root [label=\"\",shape=point,style=invis]\n";
		os << "t [label=<<font point-size=\"20\">1</font>>,shape=box,tooltip=\"1\",width=0.3,height=0.3]\n";
		
		auto toplabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u;
		auto mag = thicknessFromMagnitude(e.w);
		
		os << "root->";
		if (dd::Package::isTerminal(e)) {
			os << "t";
		} else {
			os << toplabel;
		}
		os << "[penwidth=\"" << mag <<  "\",tooltip=\"" << e.w << "\"";
		if(colored) {
			os << ",color=\"#" << colorFromPhase(e.w) << "\"";
		} else if (!CN::equalsOne(e.w)) {
			os << ",style=dashed";
		}
		if (edgeLabels) {
		    os << ",label=<<font point-size=\"8\">&nbsp;" << e.w << "</font>>";
		}

		os << "]\n";
		return os;
	}

	std::ostream& matrixNodeMatrixAndXlabel(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[label=<";
		os << R"(<font point-size="6"><table border="1" cellspacing="0" style="rounded"><tr>)";
		os << R"(<td port="0" tooltip=")" << e.p->e[0].w << R"(" href="javascript:;" sides="RB">)" << (CN::equalsZero(e.p->e[0].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>") << "</td>";
		os << R"(<td port="1" tooltip=")" << e.p->e[1].w << R"(" href="javascript:;" sides="LB">)" << (CN::equalsZero(e.p->e[1].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>") << "</td></tr><tr>";
		os << R"(<td port="2" tooltip=")" << e.p->e[2].w << R"(" href="javascript:;" sides="RT">)" << (CN::equalsZero(e.p->e[2].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>") << "</td>";
		os << R"(<td port="3" tooltip=")" << e.p->e[3].w << R"(" href="javascript:;" sides="LT">)" << (CN::equalsZero(e.p->e[3].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>") << "</td>";
		os << "</tr></table></font>>,tooltip=\"q" << e.p->v << "\"" << R"(,xlabel=<<font point-size="8">q<sub><font point-size="6">)" << e.p->v << "</font></sub></font>>]\n";
		return os;
	}

	std::ostream& matrixNodeMiddleVar(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[label=<";
		os << R"(<font point-size="10"><table border="1" cellspacing="0" cellpadding="2" style="rounded">)";
		os << R"(<tr><td colspan="2" rowspan="2" port="0" href="javascript:;" border="0" tooltip=")" << e.p->e[0].w << "\">" << (CN::equalsZero(e.p->e[0].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>")
		<< R"(</td><td sides="R"></td><td sides="L"></td>)"
		<< R"(<td colspan="2" rowspan="2" port="1" href="javascript:;" border="0" tooltip=")" << e.p->e[1].w << "\">" << (CN::equalsZero(e.p->e[1].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>")<< R"(</td></tr>)";
		os << R"(<tr><td sides="R"></td><td sides="L"></td></tr>)";
		os << R"(<tr><td colspan="2" sides="B"></td><td colspan="2" rowspan="2" border="0"><font point-size="24">q<sub><font point-size="16">)" << e.p->v << R"(</font></sub></font></td><td colspan="2" sides="B"></td></tr>)";
		os << R"(<tr><td sides="T" colspan="2"></td><td sides="T" colspan="2"></td></tr>)";
		os << R"(<tr><td colspan="2" rowspan="2" port="2" href="javascript:;" border="0" tooltip=")" << e.p->e[2].w << "\">" << (CN::equalsZero(e.p->e[2].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>")
		   << R"(</td><td sides="R"></td><td sides="L"></td>)"
		   << R"(<td colspan="2" rowspan="2" port="3" href="javascript:;" border="0" tooltip=")" << e.p->e[3].w << "\">" << (CN::equalsZero(e.p->e[3].w) ? "&nbsp;0 " : "<font color=\"white\">&nbsp;0 </font>")<< "</td></tr>";
		os << R"(<tr><td sides="R"></td><td sides="L"></td></tr>)";
		os << "</table></font>>,tooltip=\"q" << e.p->v << "\"]\n";
		return os;
	}

	std::ostream& classicMatrixNode(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[shape=circle, width=0.53, fixedsize=true, label=<";
		os << R"(<font point-size="6"><table border="0" cellspacing="0" cellpadding="0">)";
		os << R"(<tr><td colspan="4"><font point-size="18">q<sub><font point-size="10">)" << e.p->v << R"(</font></sub></font></td></tr><tr>)";
		os << R"(<td port="0" tooltip=")" << e.p->e[0].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[0].w)) {
			os << R"(<font point-size="8">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}
		os << "</td>";
		os << "<td></td><td></td>";
		os << R"(<td port="3" tooltip=")" << e.p->e[3].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[3].w)) {
			os << R"(<font point-size="8">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}		os << "</td>";
		os << "</tr><tr><td></td>";
		os << R"(<td port="1" tooltip=")" << e.p->e[1].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[1].w)) {
			os << R"(<font point-size="8">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}		os << "</td>";
		os << R"(<td port="2" tooltip=")" << e.p->e[2].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[2].w)) {
			os << R"(<font point-size="8">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}		os << "</td>";
		os << "<td></td></tr></table></font>>,tooltip=\"q" << e.p->v << "\"]\n";
		return os;
	}

	std::ostream& vectorNode(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[label=<";
		os << R"(<font point-size="8"><table border="1" cellspacing="0" cellpadding="0" style="rounded">)";
		os << R"(<tr><td colspan="2" border="0" cellpadding="1"><font point-size="20">q<sub><font point-size="12">)" << e.p->v << R"(</font></sub></font></td></tr><tr>)";
		os << R"(<td height="6" width="14" port="0" tooltip=")" << e.p->e[0].w << R"(" href="javascript:;" sides="RT">)" << (CN::equalsZero(e.p->e[0].w) ? "&nbsp;0 " : R"(<font color="white">&nbsp;0 </font>)") << "</td>";
		os << R"(<td height="6" width="14" port="2" tooltip=")" << e.p->e[2].w << R"(" href="javascript:;" sides="LT">)" << (CN::equalsZero(e.p->e[2].w) ? "&nbsp;0 " : R"(<font color="white">&nbsp;0 </font>)") << "</td>";
		os << "</tr></table></font>>,tooltip=\"q" << e.p->v << "\"]\n";
		return os;
	}

	std::ostream& vectorNodeVectorLook(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[label=<";
		os << R"(<font point-size="10"><table border="1" cellspacing="0" cellpadding="2" style="rounded">)";
		os << R"(<tr><td rowspan="2" sides="R" cellpadding="2"><font point-size="18">q<sub><font point-size="12">)" << e.p->v << "</font></sub></font></td>";
		os << R"(<td port="0" tooltip=")" << e.p->e[0].w << R"(" href="javascript:;" sides="LB">)" << (CN::equalsZero(e.p->e[0].w) ? "&nbsp;0 " : R"(<font color="white">&nbsp;0 </font>)") << "</td></tr><tr>";
		os << R"(<td port="2" tooltip=")" << e.p->e[2].w << R"(" href="javascript:;" sides="LT">)" << (CN::equalsZero(e.p->e[2].w) ? "&nbsp;0 " : R"(<font color="white">&nbsp;0 </font>)") << "</td>";
		os << "</tr></table></font>>,tooltip=\"q" << e.p->v << "\"]\n";
		return os;
	}

	std::ostream& classicVectorNode(const Edge& e, std::ostream& os) {
		auto nodelabel = ((uintptr_t)e.p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
		os << nodelabel << "[shape=circle, width=0.46, fixedsize=true, label=<";
		os << R"(<font point-size="6"><table border="0" cellspacing="0" cellpadding="0">)";
		os << R"(<tr><td colspan="2"><font point-size="18">q<sub><font point-size="10">)" << e.p->v << R"(</font></sub></font></td></tr><tr>)";
		os << R"(<td port="0" tooltip=")" << e.p->e[0].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[0].w)) {
			os << R"(<font point-size="10">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}
		os << "</td>";
		os << R"(<td port="2" tooltip=")" << e.p->e[2].w << R"(" href="javascript:;">)";
		if (CN::equalsZero(e.p->e[2].w)) {
			os << R"(<font point-size="10">&nbsp;0 </font>)";
		} else {
			os << R"(<font color="white">&nbsp;0 </font>)";
		}
		os << "</td>";
		os << "</tr></table></font>>,tooltip=\"q" << e.p->v << "\"]\n";
		return os;
	}

	std::ostream& matrixEdge(const Edge& from, const Edge& to, short idx, std::ostream& os, bool edgeLabels, bool classic, bool colored) {
		auto fromlabel = ((uintptr_t)from.p & 0x001fffffu) >> 1u;
		auto tolabel = ((uintptr_t)to.p & 0x001fffffu) >> 1u;

		os << fromlabel << ":" << idx << ":";
		if (classic) {
			if (idx == 0) os << "sw";
			else if (idx == 1 || idx == 2) os << "s";
			else os << "se";
		} else {
			if (idx == 0) os << "sw";
			else if (idx == 1) os << "se";
			else os << 's';
		}
		os << "->";
		if (Package::isTerminal(to)) {
			os << "t";
		} else {
			os << tolabel;
			if (!classic)
				os << ":n";
		}

		auto mag = thicknessFromMagnitude(to.w);
		os << "[penwidth=\"" << mag << "\",tooltip=\"" << to.w << "\"";
		if(colored) {
			os << "color=\"#" << colorFromPhase(to.w) << "\"";
		} else if (!CN::equalsOne(to.w)) {
			os << ",style=dashed";
		}
		if (edgeLabels) {
			os << ",label=<<font point-size=\"8\">&nbsp;" << to.w << "</font>>";
		}
		os << "]\n";

		return os;
	}

	std::ostream& vectorEdge(const Edge& from, const Edge& to, short idx, std::ostream& os, bool edgeLabels, bool classic, bool colored) {
		auto fromlabel = ((uintptr_t)from.p & 0x001fffffu) >> 1u;
		auto tolabel = ((uintptr_t)to.p & 0x001fffffu) >> 1u;

		os << fromlabel << ":" << idx << ":";
		os << (idx == 0 ? "sw" : "se") << "->";
		if (Package::isTerminal(to)) {
			os << "t";
		} else {
			os << tolabel;
		}

		auto mag = thicknessFromMagnitude(to.w);
		os << "[penwidth=\"" << mag << "\",tooltip=\"" << to.w << "\"";
		if(colored) {
			os <<  ",color=\"#" << colorFromPhase(to.w) << "\"";
		} else if (!CN::equalsOne(to.w)) {
			os << ",style=dashed";
		}
		if (edgeLabels) {
			os << ",label=<<font point-size=\"8\">&nbsp;" << to.w << "</font>>";
		}
		os << "]\n";

		return os;
	}

	void toDot(const Edge& e, std::ostream& os, bool isVector, bool colored, bool edgeLabels, bool classic) {
		std::ostringstream oss{};
		// header, root and terminal declaration

		header(e, oss, edgeLabels, colored);

		std::unordered_set<NodePtr> nodes{};
		auto priocmp = [] (const dd::Edge* left, const dd::Edge* right) { return left->p->v < right->p->v; };
		std::priority_queue<const dd::Edge*, std::vector<const dd::Edge*>, decltype(priocmp)> q(priocmp);
		q.push(&e);

		// bfs until finished
		while (!q.empty()) {
			auto node = q.top();
			q.pop();

			// base case
			if (Package::isTerminal(*node))
				continue;

			// check if node has already been processed
			auto ret = nodes.emplace(node->p);
			if (!ret.second) continue;

			// node definition as HTML-like label (href="javascript:;" is used as workaround to make tooltips work)
			if (isVector) {
				if (classic)
					classicVectorNode(*node, oss);
				else
					vectorNode(*node, oss);
			} else {
				if (classic)
					classicMatrixNode(*node, oss);
				else
					matrixNodeMiddleVar(*node, oss);
			}

			// iterate over edges in reverse to guarantee correct proceossing order
			for (short i=dd::NEDGE-1; i >= 0; --i) {
				if (isVector && i%2 != 0)
					continue;

				auto& edge = node->p->e[i];
				if (CN::equalsZero(edge.w)) {
					if (classic) {
						// potentially add zero stubs here
//						auto nodelabel = ((uintptr_t)node->p & 0x001fffffu) >> 1u; // this allows for 2^20 (roughly 1e6) unique nodes
//						oss << nodelabel << "0" << i << "[label=<<font point-size=\"6\">0</font>>]\n";
//						oss << nodelabel << ":" << i << "->" << nodelabel << "0" << i << ":n\n";
						continue;
					} else {
						continue;
					}
				}

				// non-zero edge to be included
				q.push(&edge);

				if (isVector) {
					vectorEdge(*node, edge, i, oss, edgeLabels, classic, colored);
				} else {
					matrixEdge(*node, edge, i, oss, edgeLabels, classic, colored);
				}
			}
		}
		oss << "}\n";

		os << oss.str() << std::flush;
	}

	void export2Dot(Edge basic, const std::string& outputFilename, bool isVector, bool colored, bool edgeLabels, bool classic, bool show) {
		std::ofstream init(outputFilename);
		toDot(basic, init, isVector, colored, edgeLabels, classic);
		init.close();

		if (show) {
			std::ostringstream oss;
			oss << "dot -Tsvg " << outputFilename << " -o " << outputFilename.substr(0, outputFilename.find_last_of('.')) << ".svg";
			auto str = oss.str(); // required to avoid immediate deallocation of temporary
			static_cast<void>(!std::system(str.c_str())); // cast and ! just to suppress the unused result warning
		}
	}

	RGB hlsToRGB(const fp& h, const fp& l, const fp& s) {
		if (s == 0.0) {
			return {l, l, l};
		}
		fp m2;
		if (l <= 0.5) {
			m2 = l * (1+s);
		} else {
			m2 = l+s-(l*s);
		}
		auto m1 = 2*l - m2;

		auto v = [] (const fp& m1, const fp& m2, fp hue) -> fp {
			while (hue < 0) hue += 1.0;
			while (hue > 1) hue -= 1.0;
			if (hue < 1./6)
				return m1 + (m2-m1)*hue*6.0;
			if (hue < 0.5)
				return m2;
			if (hue < 2./3)
				return m1 + (m2-m1)*(2./3-hue)*6.0;
			return m1;
		};

		return {v(m1, m2, h+1./3), v(m1, m2, h), v(m1, m2, h-1./3)};
	}

	RGB colorFromPhase(const Complex& a) {
		auto phase = CN::arg(a);
		auto twopi = 2*PI;
		phase = (phase) / twopi;
		return hlsToRGB(phase, 0.5, 0.5);
	}

	fp thicknessFromMagnitude (const Complex& a) {
		return 3.0*std::max(CN::mag(a), 0.10);
	}
	
	// 0 1 (1 0.7071067811865476) () () (2 0.7071067811865476+0.8i)
	void serialize(Edge basic, const std::string& outputFilename, bool isVector) {
		std::ofstream init(outputFilename);
		std::ostringstream oss{};

		serialize(basic, oss, isVector);

		init << oss.str() << std::flush;
		init.close();
	}

	void serialize(Edge basic, std::ostream& oss, bool isVector) {
		int next_index = 0;		
		std::unordered_map<NodePtr, int> node_index{};	
		
		oss << SERIALIZATION_VERSION << "\n";
		oss << CN::toString(basic.w, false, 16) << "\n";
 
		/* BFS
		std::unordered_set<NodePtr> nodes{};
		auto priocmp = [] (const dd::Edge* left, const dd::Edge* right) { return left->p->v < right->p->v; };
		std::priority_queue<const dd::Edge*, std::vector<const dd::Edge*>, decltype(priocmp)> q(priocmp);
		
		node_index[basic.p] = next_index++;	
		q.push(&basic);

		// bfs until finished
		while (!q.empty()) {
			auto node = q.top();
			q.pop();

			// base case
			if (Package::isTerminal(*node))
				continue;

			// check if node has already been processed
			auto ret = nodes.emplace(node->p);
			if (!ret.second) continue;


			oss  << node_index[node->p] << " " << node->p->v;

			// iterate over edges in reverse to guarantee correct proceossing order
			for (short i = 0; i < dd::NEDGE; ++i) {
				oss << " (";
				if (!isVector || i % 2 == 0) {
					auto& edge = node->p->e[i];
					if (!CN::equalsZero(edge.w)) {
						int current_idx;
						if(Package::isTerminal(edge)) {
							current_idx = -1;
						} else {
							const auto& nodeit = node_index.find(edge.p);
							if(nodeit != node_index.end()) {
								current_idx = nodeit->second;
							} else {
								node_index[edge.p] = current_idx = next_index;
								next_index++;
							}
						}

						q.push(&edge);
						oss << current_idx << " " << CN::toString(edge.w, false, 16);
					}
				}


				oss << ")";
			}
			oss << "\n";
		}
		*/
		
		// POST ORDER TRAVERSAL USING ONE STACK   https://www.geeksforgeeks.org/iterative-postorder-traversal-using-stack/
		
		std::stack<Edge*> stack;

		Edge *node = &basic;
		if(node != nullptr && !Package::isTerminal(*node)) {
			do {
				while(node != nullptr && !Package::isTerminal(*node)) {
					for (short i=dd::NEDGE-1; i > 0; --i) {
						if (isVector && i % 2 != 0) {
							continue;
						}
						auto& edge = node->p->e[i];
						if (Package::isTerminal(edge)) {
							continue;
						}
						if (CN::equalsZero(edge.w)) {
							continue;
						}
						if(node_index.find(edge.p) != node_index.end()) {
							continue;
						}

						// non-zero edge to be included
						stack.push(&edge);
					}
					stack.push(node);
					node = &node->p->e[0];
				}
				node = stack.top();
				stack.pop();
				
				bool hasChild = false;
				for (short i = 1; i < dd::NEDGE && !hasChild; ++i) {
					if (isVector && i % 2 != 0) {
						continue;
					}
					auto& edge = node->p->e[i];
					if (CN::equalsZero(edge.w)) {
						continue;
					}
					if(node_index.find(edge.p) != node_index.end()) {
						continue;
					}
					hasChild = edge.p == stack.top()->p;
				}

				if(hasChild) {
					dd::Edge* temp = stack.top();
					stack.pop();
					stack.push(node);
					node = temp;
				} else {
					if(node_index.find(node->p) != node_index.end()) {
						node = nullptr;
						continue;
					}
					node_index[node->p] = next_index;
					next_index++;
					oss << node_index[node->p] << " " << node->p->v;

					// iterate over edges in reverse to guarantee correct processing order
					for (short i = 0; i < dd::NEDGE; ++i) {
						oss << " (";
						if (isVector || i % 2 == 0) {
							auto& edge = node->p->e[i];
							if (!CN::equalsZero(edge.w)) {
								int edge_idx = Package::isTerminal(edge) ? -1 : node_index[edge.p];			
								oss << edge_idx << " " << CN::toString(edge.w, false, 16);
							}
						}
						oss << ")";
					}
					oss << "\n";
					node = nullptr;
				}
			} while (!stack.empty());
		}
		/* POST ORDER TRAVERSAL USING TWO STACKS   https://www.geeksforgeeks.org/iterative-postorder-traversal/

		std::stack<Edge*> stack;
		std::unordered_set<NodePtr> nodes;
		stack.push(&basic);
		std::stack<Edge*> reverse_stack;
		// iterative post order traversal
		while (!stack.empty()) {
			auto node = stack.top();
			stack.pop();
			
			// base case
			if (Package::isTerminal(*node))
				continue;

			// check if node has already been processed
			auto ret = nodes.emplace(node->p);
			if (!ret.second) continue;

			reverse_stack.push(node);
			for (short i=dd::NEDGE-1; i >= 0; --i) {
				if (isVector && i % 2 != 0)
					continue;

				auto& edge = node->p->e[i];
				if (CN::equalsZero(edge.w)) {
					continue;
				}

				// non-zero edge to be included
				stack.push(&edge);
			}
		}

		while(!reverse_stack.empty()) {
			auto node = reverse_stack.top();
			reverse_stack.pop();

			int current_idx;
			const auto& nodeit = node_index.find(node->p);
			if(nodeit != node_index.end()) {
				current_idx = nodeit->second;
			} else {
				node_index[node->p] = current_idx = next_index;
				next_index++;
			}
			oss << "n: " << current_idx << " " << node->p->v << "\n";

			// iterate over edges in reverse to guarantee correct proceossing order
			for (short i=dd::NEDGE-1; i >= 0; --i) {
				if (isVector && i % 2 != 0)
					continue;
				auto& edge = node->p->e[i];
				if (CN::equalsZero(edge.w)) {
					continue;
				}
				int edge_idx = Package::isTerminal(edge) ? -1 : node_index[edge.p];

				auto real = ComplexNumbers::val(edge.w.r);
				auto imag = ComplexNumbers::val(edge.w.i);
				oss << "[" << i << "]: " << edge_idx << " " << real << " " << imag << "\n";
			}
		}		
		*/
	}

	dd::Edge create_deserialized_node(std::unique_ptr<dd::Package>& dd, int index, int v, std::array<int, dd::NEDGE>& edge_idx, 
									  std::array<dd::ComplexValue, dd::NEDGE>& edge_weight, std::unordered_map<int, NodePtr>& nodes) {
		if(index == -1) {
			return dd::Package::DDzero;
		}

		std::array<dd::Edge, dd::NEDGE> edges{};
		for(int i = 0; i < dd::NEDGE; i++) {
			if(edge_idx[i] == -2) {
				edges[i] = dd::Package::DDzero;
			} else {
				edges[i].p = edge_idx[i] == -1 ? dd::Package::DDone.p : nodes[edge_idx[i]];
				edges[i].w = dd->cn.lookup(edge_weight[i]);
			}
		}

		dd::Edge newedge = dd->makeNonterminal(v, edges);
		nodes[index] = newedge.p;
		
		// reset
		edge_idx.fill(-2);

		return newedge;
	}

	ComplexValue toComplexValue(const std::string& real_str, std::string imag_str) {
		fp real = real_str.empty() ? 0 : std::stod(real_str); 
		
		imag_str.erase(remove(imag_str.begin(), imag_str.end(), ' '), imag_str.end());
		imag_str.erase(remove(imag_str.begin(), imag_str.end(), 'i'), imag_str.end());
		if(imag_str == "+" || imag_str == "-") imag_str = imag_str + "1";
		fp imag = imag_str.empty() ? 0 : std::stod(imag_str); 
		return dd::ComplexValue{real, imag};
	}

	dd::Edge deserialize(std::unique_ptr<dd::Package>& dd, const std::string& inputFilename) {
		auto ifs = std::ifstream(inputFilename);
		
		if(!ifs.good()) {
			std::cerr << "Wrong Version of serialization file" << std::endl;
			exit(1);
		}

		return deserialize(dd, ifs);
	}

	dd::Edge deserialize(std::unique_ptr<dd::Package>& dd, std::istream& ifs) {
		std::string version;
		std::getline(ifs, version);
		// ifs >> version;
		if(strcmp(version.c_str(), SERIALIZATION_VERSION) != 0) {
			std::cerr << "Wrong Version of serialization file version: " << version << std::endl;
			exit(1);
		}
		
		dd::Edge result = dd::Package::DDzero;
		std::unordered_map<int, NodePtr> nodes;
		
		std::string line;		
		std::string complex_real_regex = R"(([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![ \d\.]*(?:[eE][+-])?\d*[iI]))?)";
    	std::string complex_imag_regex = R"(( ?[+-]? ?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?[iI])?)";
    	std::string edge_regex = " \\(((-?\\d+) (" + complex_real_regex + complex_imag_regex + "))?\\)";
		std::regex complex_weight_regex (complex_real_regex + complex_imag_regex);
    	std::regex line_regex ("(\\d+) (\\d+)(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ")(?:" + edge_regex + ") *(?:#.*)?");
    	// std::regex e ("(\\d+) (\\d+)(?:" + edge_regex + "){4} *#.*"); // TODO {4} overwrites groups
		std::smatch m;
		
		ComplexValue rootweight{};
		int node_index;
		int v;
		std::array<int, dd::NEDGE> edge_indices{};
		edge_indices.fill(-2);
		std::array<dd::ComplexValue, dd::NEDGE> edge_weights{};

		if(std::getline(ifs, line)) {
			if(!std::regex_match(line, m, complex_weight_regex)) {
				std::cerr << "Regex did not match second line: " << line << std::endl;
				exit(1);
			}
			rootweight = toComplexValue(m.str(1), m.str(2));
		}
		// std::cout << "rootweight = " << rootweight << std::endl;

		while (std::getline(ifs, line)) {
			if (line.empty() || line.size() == 1) continue;

    		if(!std::regex_match(line, m, line_regex)) {
				std::cerr << "Regex did not match line: " << line << std::endl;
				exit(1);
			}

			// match 1: node_idx 
			// match 2: qubit_idx 

			// repeats for every edge
			// match 3: edge content
			// match 4: edge_target_idx 
			// match 5: real + imag (without i)
			// match 6: real
			// match 7: imag (without i)
			node_index = std::stoi(m.str(1));
			v          = std::stoi(m.str(2));
			 
			// std::cout << "nidx: " << node_index << " v: " << v << std::endl;
			
			for(int edge_idx = 3, i = 0; i < dd::NEDGE; i++, edge_idx += 5) {
				if(m.str(edge_idx).empty()) {
					// std::cout << "index " << i << " is empty " << std::endl;
					continue;
				}

				edge_indices[i] = std::stoi(m.str(edge_idx + 1));
				edge_weights[i] = toComplexValue(m.str(edge_idx + 3), m.str(edge_idx + 4)); 
			}

			result = create_deserialized_node(dd, node_index, v, edge_indices, edge_weights, nodes);
		}

		Complex w = dd->cn.getCachedComplex(rootweight.r, rootweight.i);
		CN::mul(w, result.w, w);
		result.w = dd->cn.lookup(w);
		dd->cn.releaseCached(w);

        return result;
		/*
		auto ifs = std::ifstream(inputFilename);
		if(!ifs.good()) {
			std::cerr << "Wrong Version of serialization file" << std::endl;
			exit(1);
		}
	
		std::string version;
		ifs >> version;
		if(strcmp(version.c_str(), SERIALIZATION_VERSION) != 0) {
			std::cerr << "Wrong Version of serialization file" << std::endl;
			exit(1);
		}
		
		std::string line;
		std::string prefix;
		
		int index = -1, v = -1;
		std::array<int, dd::NEDGE> edge_idx;
		edge_idx.fill(-2);
		std::array<dd::ComplexValue, dd::NEDGE> edge_weight;
		fp real, imag;
		std::unordered_map<int, NodePtr> nodes;

		while (std::getline(ifs, line)) {
			if (line.empty() || line.size() == 1) continue;
			
			std::stringstream ss(line);
			ss >> prefix;
			if(prefix == "n:") {
				create_deserialized_node(dd, index, v, edge_idx, edge_weight, nodes);
				ss >> index;
				ss >> v;
			} else {
				int i = static_cast<int>(prefix[1] - '0');
				ss >> edge_idx[i];
				ss >> real;
				ss >> imag;
				edge_weight[i] = {real, imag};
			}
		}
		std::cout << "last create" << std::endl;
		return create_deserialized_node(dd, index, v, edge_idx, edge_weight, nodes);
		*/
	}
}
