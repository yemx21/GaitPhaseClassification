#include "FeatureExtraction.h"
#include <iostream>
#include "P1D.hpp"
#include "FileStream.h"

template <typename It>
auto Median(It begin, It end)
{
	const auto size = std::distance(begin, end);
	std::nth_element(begin, begin + size / 2, end);
	return *std::next(begin, size / 2);
}

struct PairHashFun
{
	size_t operator()(const std::pair<int, int>& x) const throw()
	{
		return std::hash<int>()(x.first) ^ std::hash<int>()(x.second);
	}
};

class PairEqualFun
{
public:
	bool operator() (const std::pair<int, int>& x1, const std::pair<int, int>& x2) const
	{
		return (x1.first == x2.second && x1.second == x2.first) || (x1.first == x2.first && x1.second == x2.second);
	}
};

struct FeatureQuality
{
public:
	sizetype Major;
	double Gini;
	double Importance;
	std::map<sizetype, double> Class;
	

	FeatureQuality() : Gini(0), Importance(1.0)
	{

	}

	FeatureQuality(double gini, double impo, const std::map<sizetype, double>& cls) : Gini(gini), Importance(impo), Class(cls)
	{
		Major = -1;
		double maxp = 0.0;
		for (const auto& c : cls)
		{
			if (c.second > maxp)
			{
				maxp = c.second;
				Major = c.first;
			}
		}
	}
};

void FeatureExtractor::ExtractFeaturePairs(SequenceDataSet* seqdata, const std::vector<sizetype>& channels, const std::vector<sizetype>& labels, std::vector<std::vector<std::pair<int, int>>>& finalpairs, int desiredfeatcount, int neigbor)
{
	p1d::Persistence1D p;

	int channel_count = (int)channels.size();

	std::map<int, int> chmap;
	for (int ci = 0; ci < channel_count; ci++)
	{
		chmap[channels[ci]] = ci;
	}

	std::vector<std::unordered_map<std::pair<int, int>, long long, PairHashFun, PairEqualFun>> pairs;
	pairs.resize(channels.size());

	auto samples = seqdata->getTemplates();

	std::map<sizetype, long double> temcls;
	for (sizetype l : labels) temcls[l] = 0.0;

	const auto& sqam = samples.begin()->second;
	{
		for (auto& s : sqam)
		{
			for (auto l : s->Labels)
			{
				temcls[l] += 1.0;
			}
		}
	}

	long double stemcls = 0.0;
	for (const auto& s : temcls)
	{
		stemcls += s.second;
	}

	std::map<sizetype, double> wcls;
	for (const auto& s : temcls)
	{
		wcls[s.first] = (double)(s.second / stemcls);
	}

	for (auto& sam : samples)
	{
		int c = sam.first;
		auto& cpairs = pairs[chmap[c]];
		for(auto& s: sam.second)
		{
			const auto& sdata = s->Values;

			int datalen = sdata.size();
			int halflen = (int)(datalen / 2);

			//copy left half, right half
			std::vector<float> chdata = sdata;
			chdata.insert(chdata.begin(), sdata.begin() + halflen, sdata.end());
			chdata.insert(chdata.end(), sdata.begin(), sdata.begin() + halflen);
			
			p.RunPersistence(chdata, -1);

			std::vector<int> tmins;
			std::vector<int> tmaxs;
			p.GetExtremaIndices(tmins, tmaxs);

			std::vector<p1d::TPairedExtrema> textremas;
			p.GetPairedExtrema(textremas);

			std::set<int> mins;
			std::set<int> maxs;

			for (int m : tmins)
			{
				if (m < halflen)
					mins.insert(halflen + m);
				else if(m < halflen + datalen)
					mins.insert(m - halflen);
				else
					mins.insert(m- halflen - datalen);
			}
			for (int m : tmaxs)
			{
				if (m < halflen)
					maxs.insert(halflen + m);
				else if (m < halflen + datalen)
					maxs.insert(m - halflen);
				else
					maxs.insert(m - halflen - datalen);
			}

			std::vector<p1d::TPairedExtrema> extremas;

			for (const auto& te : textremas)
			{
				p1d::TPairedExtrema nte;
				if (te.MinIndex < halflen)
					nte.MinIndex = halflen + te.MinIndex;
				else if (te.MinIndex < halflen + datalen)
					nte.MinIndex = te.MinIndex - halflen;
				else
					nte.MinIndex = te.MinIndex - halflen - datalen;

				if (te.MaxIndex < halflen)
					nte.MaxIndex = halflen + te.MaxIndex;
				else if (te.MaxIndex < halflen + datalen)
					nte.MaxIndex = te.MaxIndex - halflen;
				else
					nte.MaxIndex = te.MaxIndex - halflen - datalen;

				nte.Persistence = te.Persistence;
				extremas.push_back(nte);
			}

			std::vector<float> betas;
			std::vector<int> betaindices;
			for (int minidx : mins)
			{
				for (const auto& ex : extremas)
				{
					if (ex.MinIndex == minidx || ex.MaxIndex == minidx)
					{
						betas.push_back(ex.Persistence);
						betaindices.push_back(minidx);
						break;
					}
				}
			}
			for (int maxidx : maxs)
			{
				for (const auto& ex : extremas)
				{
					if (ex.MinIndex == maxidx || ex.MaxIndex == maxidx)
					{
						betas.push_back(ex.Persistence);
						betaindices.push_back(maxidx);
						break;
					}
				}
			}

			int extrenum = betas.size();
			if (extrenum)
			{
				float beta_thres = 0.0;
				if (extrenum < 20)// sqrt(time_step))
				{
					beta_thres = Median(betas.begin(), betas.end());
				}

				//std::cout << "beta_thres=" << beta_thres << std::endl;

				auto iter1 = betas.begin();
				auto iter2 = betaindices.begin();
				for (; iter1 != betas.end();)
				{
					if (*iter1 < beta_thres)
					{
						iter1 = betas.erase(iter1);
						iter2 = betaindices.erase(iter2);
					}
					else
					{
						iter1++;
						iter2++;
					}
				}

				int betanum = betas.size();

				for (int i = 0; i < betanum; i++)
				{
					for (int j = i; j < betanum; j++)
					{
						int idx1 = betaindices[i];
						int idx2 = betaindices[j];

						if (abs(idx1 - idx2) > neigbor)
						{
							std::pair<int, int> key;
							key.first = idx1;
							key.second = idx2;

							auto p = cpairs.find(key);
							if (p == cpairs.end())
								cpairs.insert(std::make_pair(key, 0ll));
						}
					}
				}

			}
		}

	}

	std::vector<std::unordered_map<std::pair<int, int>, std::map<sizetype, std::pair<double, double>>, PairHashFun, PairEqualFun>> fpairs;
	fpairs.resize(channels.size());

	for (auto& sam : samples)
	{
		int c = sam.first;
		auto& cp = pairs[chmap[c]];
		auto& fp = fpairs[chmap[c]];
		for (auto& s : sam.second)
		{
			const auto& sdata = s->Values;

			int datalen = sdata.size();
			int halflen = (int)(datalen / 2);

			//copy left half, right half
			std::vector<float> chdata = sdata;
			chdata.insert(chdata.begin(), sdata.begin() + halflen, sdata.end());
			chdata.insert(chdata.end(), sdata.begin(), sdata.begin() + halflen);

			for (int i = 0; i < datalen; i++)
			{
				for (const auto& f : cp)
				{
					int f1 = f.first.first;
					int f2 = f.first.second;
					if (f1 > f2)
					{
						f1 = f.first.second;
						f2 = f.first.first;
					}
					int U = 0;
					int V = 0;
					if (i < f1)
					{
						U = f2 - datalen - i;
						V = f2 - i;
					}
					else if (i > f2)
					{
						U = i - f2;
						V = datalen - i + f1;
					}
					else
					{
						U = i - f1;
						V = f2 - i;
					}
					
					std::pair<int, int> key;
					key.first = U;
					key.second = V;

					if (fp.find(key) == fp.end())
					{
						std::map<sizetype, std::pair<double, double>> emptyprobs;
						for (auto lab : labels)
						{
							emptyprobs[lab].first = 0;
							emptyprobs[lab].second = 0;
						}

						fp.insert(std::make_pair(key, emptyprobs));
					}
					auto & kprobs = fp[key];
					for (auto sp : s->Probabilities[i])
					{
						//kprobs[sp.first] += sp.second;// *s->LocalWeight * s->GlobalWeight;
						kprobs[sp.first].first += sp.second *s->LocalWeight * s->GlobalWeight;
						kprobs[sp.first].second += s->LocalWeight * s->GlobalWeight;
					}
				}
			}
		}
	}

	std::vector<std::unordered_map<std::pair<int, int>, FeatureQuality, PairHashFun, PairEqualFun>> qpairs;
	qpairs.resize(channels.size());

	for (int c=0; c<channel_count; c++)
	{
		const auto& fp = fpairs[c];
		auto& qp = qpairs[c];

		for (const auto& f : fp)
		{
			qp.insert(std::make_pair(f.first, FeatureQuality()));

			double np = 0.0;
			for (const auto& ff : f.second)
			{
				np += ff.second.first;
			}
			double np2 = np*np;

			double impo = 0.0;

			std::map<sizetype, double> clsp;
			for (const auto& ff : f.second)
			{
				clsp[ff.first] = ff.second.first / np;
				impo += ff.second.second;
			}
			
			double sp = 0.0;
			for (const auto& ff : f.second)
			{
				sp += (ff.second.first*ff.second.first)/np2;
			}

			qp[f.first] = FeatureQuality(sp, impo, clsp);
		}
	}

	std::vector<std::vector<std::pair<std::pair<int, int>, FeatureQuality>>> opairs;
	std::vector<double> mpairs;

	double mpairsum = 0.0;
	opairs.resize(channel_count);

	for (int c = 0; c < channel_count; c++)
	{
		const auto& q = qpairs[c];
		auto& copairs = opairs[c];
		std::vector<double> occurs;
		double occurnum = 0.0;
		for (const auto& cp : q)
		{
			occurnum += cp.second.Gini;
		}
		for (const auto& cp : q)
		{
			std::pair<std::pair<int, int>, FeatureQuality> cfp;
			cfp.first = cp.first;
			cfp.second = FeatureQuality((double)cp.second.Gini, cp.second.Importance, cp.second.Class);
			//cfp.second = FeatureQuality((double)cp.second.Gini / occurnum, cp.second.Class);
			occurs.push_back(cfp.second.Gini);
			copairs.push_back(cfp);
		}
		mpairs.push_back(Median(occurs.begin(), occurs.end()));
		mpairsum += mpairs.back();
	}

	std::map<sizetype, std::vector<double>> cdb;
	for (sizetype l : labels) cdb[l].resize(channel_count);

	for (int c = 0; c < channel_count; c++)
	{
		const auto& copairs = opairs[c];

		for (const auto& cop : copairs)
		{
			cdb[cop.second.Major][c] += cop.second.Gini;
		}
	}

	std::map<sizetype, double> ctb;
	double sctb = 0.0;
	for (const auto& c : cdb)
	{
		double sc = 0.0;
		for (const auto& cc : c.second)
		{
			sc += cc;
		}
		ctb[c.first] = sc;
		sctb += sc;
	}

	finalpairs.resize(channel_count);
	if (desiredfeatcount == -1)
	{
		for (int c = 0; c < channel_count; c++)
		{
			auto& copairs = opairs[c];
			auto& fpairs = finalpairs[c];

			std::sort(std::begin(copairs), std::end(copairs), [](const auto& lhs, const auto& rhs) {
				return lhs.second.Gini > rhs.second.Gini;
			});

			for (int n = 0; n < copairs.size(); n++)
			{
				fpairs.push_back(copairs[n].first);
			}
		}
		return;
	}

	std::map<sizetype, int> selclsnums;
	for (const auto& c : ctb)
	{
		selclsnums[c.first] = (int)((double)desiredfeatcount * wcls[c.first]);
	}

	std::vector<std::map<sizetype, int>> selnums;
	selnums.resize(channel_count);

	for (int c = 0; c < channel_count; c++)
	{
		auto& copairs = opairs[c];
		auto& fpairs = finalpairs[c];

		int cfeatcount = (int)((double)desiredfeatcount * mpairs[c] / mpairsum);
		std::sort(std::begin(copairs), std::end(copairs), [](const auto& lhs, const auto& rhs) {
			return lhs.second.Gini > rhs.second.Gini;
		});

		for (int n = 0; n < cfeatcount && n < copairs.size(); n++)
		{
			fpairs.push_back(copairs[n].first);
		}
		fpairs.push_back(std::make_pair<int,int>(0,0));
	}

	/*int ch = 1;
	for (auto fp : finalpairs)
	{
		std::cout << "channel " << ch++ << std::endl;
		for (auto pp : fp)
		{
			std::cout << pp.first << "," << pp.second << std::endl;
		}
	}*/

}
