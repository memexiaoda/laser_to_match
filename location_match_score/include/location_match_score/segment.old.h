#ifndef SEGMENT_H_
#define SEGMENT_H_

#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

using namespace std;

template <typename T>
vector<vector<T> > segment(const vector<T>& input)
{
	vector<vector<T> > p_set;

	T pre_pp = input[0];

	std::cout << "input size = " << input.size() << std::endl;

	// std::cout << "first element: " << input[0] << std::endl;
	// std::cout << "pre_pp.norm() = " << pre_pp.norm() << std::endl;

	vector<T> temp;

	for(int i = 1; i < input.size(); ++i)
	{

		if(std::sqrt(input[i].SquaredDistance(pre_pp)) < 0.3)
		{
			// std::cout << "ok " << i << std::endl;
			temp.push_back(input[i]);

		}else{
			
      if(temp.size() > 10)
      {
      	p_set.push_back(temp);
      	temp.clear();
      }
      temp.clear();
		}

		pre_pp = input[i];
	}

	if(temp.size() > 10)
		p_set.push_back(temp);

	std::cout << "p_set size = " << p_set.size() << std::endl;

	// if(p_set.back().size() < 10)
	// 	p_set.back().clear();

  return p_set;
}

#endif