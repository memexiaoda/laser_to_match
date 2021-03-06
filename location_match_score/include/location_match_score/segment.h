#ifndef SEGMENT_H_
#define SEGMENT_H_

#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <yaml-cpp/yaml.h>

using namespace std;

template <typename T>
vector<vector<T> > segment(const vector<T>& input, float dis_threshold, int batch_size)
{
	vector<vector<T> > p_set;
	T pre_pp = input[0];
	// std::cout << "input size = " << input.size() << std::endl;
	vector<T> temp;
	for(int i = 1; i < input.size(); ++i)
	{
		if(std::sqrt(input[i].SquaredDistance(pre_pp)) < dis_threshold)
		{
			// std::cout << "ok " << i << std::endl;
			temp.push_back(input[i]);
		}else{
			if(temp.size() > batch_size)
			{
				p_set.push_back(temp);
				temp.clear();
			}
			temp.clear();
	  }
		pre_pp = input[i];
    }

	if(temp.size() > batch_size)
		p_set.push_back(temp);

	// std::cout << "p_set size = " << p_set.size() << std::endl;
	// if(p_set.back().size() < 10)
	// 	p_set.back().clear();
  return p_set;
}

#endif