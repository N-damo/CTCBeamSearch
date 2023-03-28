/*
** Description：beam search are revised based on fast_ctc_decode(https://github.com/nanoporetech/fast-ctc-decode) following MIT license.
** Author：linlian
** Date:
** Modify Record:
*/


#include "beam_search.hpp"
#include <random>

using namespace std;



int main()
{
	string alphabet = "NACGT";
	size_t beam_size = 5;
	float beam_cut_threshold = 0.1f;
	bool collapse_repeats = true;
	pair<string, vector<int>> r;

	constexpr int n_steps_in = 1200;
	constexpr int alphabet_size = 5;

	float posteriors[n_steps_in][alphabet_size];

	//std::random_device rd;
	std::mt19937 gen(0);
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	// Generate random posteriors
	for (int i = 0; i < n_steps_in; ++i) {
		for (int j = 0; j < alphabet_size; ++j) {
			posteriors[i][j] = dis(gen);
		}
	}

	clock_t start, end;

	start = clock();
	beam_search(posteriors[0], n_steps_in, alphabet, beam_size, beam_cut_threshold, collapse_repeats, r);
	end = clock();
	auto endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time:" << endtime * 1000000 << " us" << endl;

	if(r.first.size() != 0);
	{
		cout << r.first << endl;
		for (auto& d : r.second)
		{
			printf("%d,", d);
		}
		cout << endl;
	}
	return 0;
}
