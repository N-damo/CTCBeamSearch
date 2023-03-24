/*
** (C) Copyright 广州孔确基因科技有限公司 All Rights Reserved.
** Description：
** Author：linlian
** Date:
** Modify Record:
*/

#include <stddef.h>
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <optional>
#include <cstring>
#include <vector>
#include <cstdlib>

using namespace std;

const int ROOT_NODE = -1;






template<typename T>
void truncate_beam(vector<T>& beam, size_t beam_size)
{
	//cout << "beam size" << beam_size << endl;
	if (beam.size() > beam_size)
	{
		beam.resize(beam_size);
	}
}



struct SearchPoint
{
	int node;
	size_t state;
	float label_prob;
	float gap_prob;


	float probability()
	{
		return label_prob + gap_prob;
	}
};

template<typename T>
struct LabelNode
{
	// The index into the alphabet of this label.
	//
	// Note that blanks are not represented by a LabelNode - this is an actual label.
	size_t label;
	// The index of the parent LabelNode.
	int parent;
	// Extra data attached to the node
	T data;
};



template <typename T>
class SuffixTreeIter
{
public:
	SuffixTreeIter(const std::vector<LabelNode<T>>& nodes, int next_node)
		: nodes(nodes), next_node(next_node)
	{
	}

	std::pair<size_t, const T*> next()
	{
		if (next_node >= 0)
		{
			const LabelNode<T>& node = nodes[next_node];
			next_node = node.parent;
			return std::make_pair(node.label, &node.data);
		}
		else
		{
			return std::make_pair(0, nullptr);
		}
	}

private:
	const std::vector<LabelNode<T>>& nodes;
	int next_node;
};


template <typename T>
class SuffixTree
{

private:
	vector<LabelNode<T>> nodes;
	vector<vector<int>> children;
	vector<int> root_children;
	size_t alphabet_size;

public:

	SuffixTree(size_t alphabet_size) :alphabet_size(alphabet_size)
	{
		//this->root_children = { -1,alphabet_size };
		for (int i = 0; i < alphabet_size; i++)
		{
			this->root_children.push_back(-1);
		}

	}

	optional<size_t> label(int node)
	{
		if (node >= 0)
		{
			return this->nodes[node].label;
		}
		else
		{
			return std::nullopt;
		}
	}

	optional<int> get_child(int node, size_t label)
	{
		//print_vector(cout, root_children, ",");
		//cout <<"root children"<< endl;
		if (node == ROOT_NODE)
		{
			auto idx = this->root_children[label];
			//cout << "get child:" << idx << endl;
			if (idx >= 0)
			{
				return idx;
			}
		}
		else
		{
			assert(node >= 0);
			auto idx = this->children[node][label];
			if (idx >= 0)
			{
				return idx;
			}
		}
		return std::nullopt;
	}


	optional<int> add_node(int parent, size_t label, int data)
	{

		assert(label < this->root_children.size());
		assert(this->nodes.size() < (size_t)(INT_MAX));

		auto new_node_idx = (int)(this->nodes.size());
		if (parent == ROOT_NODE)
		{
			assert(this->root_children[label] == -1);
			this->root_children[label] = new_node_idx;
		}
		else
		{
			assert(parent >= 0);
			assert(this->children[parent][label] == -1);
			this->children[parent][label] = new_node_idx;
		}
		this->nodes.push_back(LabelNode<T>{
			label,
				parent,
				data,
		});
		this->children.push_back(std::vector<int>(alphabet_size, -1));
		return new_node_idx;
	}



	SuffixTreeIter<T> iter_from(int node) const
	{
		assert(static_cast<size_t>(node) < this->nodes.size());
		return SuffixTreeIter<T>(nodes, node);
	}

};

bool compareByNode(const SearchPoint& a, const SearchPoint& b)
{
	return a.node < b.node;
}



void beam_search(const float* network_output, const int n_steps_in,
	const string& alphabet,
	size_t beam_size,
	float beam_cut_threshold,
	bool collapse_repeats, pair<string, vector<int>>& r)
{
	vector<int> path;
	string sequence;

	size_t alphabet_size = alphabet.size();



	vector<SearchPoint> next_beam;
	vector<SearchPoint> beam;
	beam.push_back(SearchPoint{ ROOT_NODE,0,0.0f,1.0 });

	// alphabet size minus the blank label
	size_t bases = alphabet_size - 1;
	SuffixTree<int> suffix_tree(bases);

	for (int idx = 0; idx < n_steps_in; idx++)
	{
		auto pr = network_output + idx * alphabet_size;
		next_beam.clear();

		for (auto& sp : beam)
		{
			auto node = sp.node;
			auto label_prob = sp.label_prob;
			auto gap_prob = sp.gap_prob;
			auto state = sp.state;

			auto tip_label = suffix_tree.label(node);
			if (pr[0] > beam_cut_threshold)
			{
				next_beam.push_back(SearchPoint{ node,state,0.0f,(label_prob + gap_prob) * pr[0] });
			}


			for (size_t label = 0; label < bases; label++)
			{
				float pr_b = pr[label + 1];
				if (pr_b < beam_cut_threshold)
				{
					continue;
				}
				if (collapse_repeats && tip_label.has_value() && tip_label.value() == label)
				{
					next_beam.push_back(SearchPoint{
						node,
						state,
						label_prob * pr_b,
						0.0f,

						});

					auto new_node_idx = suffix_tree.get_child(node, label);
					if (!new_node_idx.has_value())
					{
						if (gap_prob > 0.0f)
						{
							new_node_idx = suffix_tree.add_node(node, label, idx);
						}
						else
						{
							new_node_idx = nullopt;
						}

					}

					if (new_node_idx.has_value())
					{

						next_beam.push_back(SearchPoint{ new_node_idx.value(),state,gap_prob * pr_b,0.0f });
					}

				}
				else
				{
					auto new_node_idx = suffix_tree.get_child(node, label);

					if (!new_node_idx.has_value())
					{
						new_node_idx = suffix_tree.add_node(node, label, idx);
					}
					next_beam.push_back(SearchPoint{
						new_node_idx.value(),
						state,
						(label_prob + gap_prob) * pr_b,
						0.0f,
						});
				}
			}

		}
		std::swap(beam, next_beam);
		const int DELETE_MARKER = INT_MIN;
		std::sort(beam.begin(), beam.end(), compareByNode);
		int last_key = DELETE_MARKER;
		int last_key_pos = 0;
		auto beam_n = beam.size();



		for (int i = 0; i < beam_n; i++)
		{
			auto& beam_item = beam[i];
			if (beam_item.node == last_key)
			{

				beam[last_key_pos].label_prob += beam_item.label_prob;
				beam[last_key_pos].gap_prob += beam_item.gap_prob;
				beam[i].node = DELETE_MARKER;
			}
			else
			{
				last_key_pos = i;
				last_key = beam_item.node;
			}
		}


		auto is_delete_marker = [&DELETE_MARKER](const SearchPoint& point) { return point.node == DELETE_MARKER; };
		auto new_end = std::remove_if(beam.begin(), beam.end(), is_delete_marker);
		beam.erase(new_end, beam.end());
		bool has_nans = false;
		std::sort(beam.begin(), beam.end(), [&](auto& a, auto& b)
			{
				auto p_a = a.probability();
				auto p_b = b.probability();
				if (std::isnan(p_a) || std::isnan(p_b))
				{
					has_nans = true;
					return false;
				}
				return p_a > p_b;
			});

		if (has_nans)
		{
			cout << "has nan" << endl;
			return;
		}

		truncate_beam(beam, beam_size);

		if (beam.size() == 0)
		{
			cout << "beam is empty" << endl;
			return;
		}

		auto top = beam[0].probability();
		for (auto& x : beam)
		{
			x.label_prob /= top;
			x.gap_prob /= top;
		}
	}



	if (beam[0].node != ROOT_NODE)
	{
		SuffixTreeIter<int> iter = suffix_tree.iter_from(beam[0].node);
		std::pair<size_t, const int*> item;
		while ((item = iter.next()).second != nullptr)
		{

			auto time = item.second;
			auto label = item.first;
			path.push_back(*(time));
			sequence += alphabet[label + 1];
		}
	}
	reverse(path.begin(), path.end());
	reverse(sequence.begin(), sequence.end());
	//cout <<"c++:" << sequence << endl;
	swap(r.first, sequence);
	swap(r.second, path);
	return;
}
