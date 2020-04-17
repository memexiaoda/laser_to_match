#ifndef _FIT_LINE_H_
#define _FIT_LINE_H_

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

namespace cv_ws{

  using namespace cv;

    template <typename T>
    double cv_getLinePara(const T& input)
    {
        //获取二维点集
        std::vector<Point> point_set;
        cv::Point point_temp;
        for( int i = 0; i < input.size(); ++i)
        {
            point_temp.x = input[i].GetX;
            point_temp.y = input[i].GetY;
            point_set.push_back(point_temp);	
        }
    
        //直线拟合 	
        //拟合结果为一个四元素的容器，比如Vec4f - (vx, vy, x0, y0)
        //其中(vx, vy) 是直线的方向向量
        //(x0, y0) 是直线上的一个点
        cv::Vec4f fitline;
        //拟合方法采用最小二乘法
        cv::fitLine(point_set, fitline, CV_DIST_L2, 0, 0.01, 0.01);
        
        //求出直线上的两个点
        double k_line = fitline[1]/fitline[0];

        return k_line;

    }
}//cv_ws



#endif