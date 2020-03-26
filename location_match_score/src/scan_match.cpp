#include <string>
#include <fstream>
#include <limits>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
// #include <boost/thread.hpp>

#include "location_match_score/scan_match.h"

namespace karto
{
ScanMatcher::~ScanMatcher()
  {
  //  delete m_pCorrelationGrid;
    // m_pCorrelationGrid = NULL;
    delete m_pSearchSpaceProbs;
    delete m_pGridLookup;

  }

  ScanMatcher* ScanMatcher::Create(Mapper* pMapper, double searchSize, double resolution,
                                   double smearDeviation, double rangeThreshold, 
                                   CorrelationGrid* pCorrelationGrid)
  {
    // invalid parameters
    if (resolution <= 0)
    {
      return NULL;
    }
    if (searchSize <= 0)
    {
      return NULL;
    }
    if (smearDeviation < 0)
    {
      return NULL;
    }
    if (rangeThreshold <= 0)
    {
      return NULL;
    }

    assert(math::DoubleEqual(math::Round(searchSize / resolution), (searchSize / resolution)));

    // calculate search space in grid coordinates
    unsigned int searchSpaceSideSize = static_cast<unsigned int>(math::Round(searchSize / resolution) + 1);
    // create search space probabilities
    Grid<double>* pSearchSpaceProbs = Grid<double>::CreateGrid(searchSpaceSideSize,
                                                                     searchSpaceSideSize, resolution);

    ScanMatcher* pScanMatcher = new ScanMatcher(pMapper);

    // std::cout<< "ok, it is now to make a grid"<<std::endl;
    pScanMatcher->m_pCorrelationGrid = pCorrelationGrid;
    std::cout << "pCorrelationGrid " << pCorrelationGrid << std::endl;

    std::cout << "pScanMatcher->m_pCorrelationGrid " << pScanMatcher->m_pCorrelationGrid << std::endl;

    uint8_t* pByte = pCorrelationGrid->GetDataPointer();

    std::cout << "pCorrelationGrid grid size = " << sizeof(pByte) << std::endl;

    pScanMatcher->m_pSearchSpaceProbs = pSearchSpaceProbs;

    std::cout<<"initialize GridLookup" <<std::endl;
    // for(int i = 0; i < 60; ++i)
    //   std::cout<< "pCorrelationGrid data i = "<< pCorrelationGrid->GetDataPointer()[i]<<std::endl;

    pScanMatcher->m_pGridLookup = new GridIndexLookup<uint8_t>(pCorrelationGrid);

    for(int i = 0; i < pCorrelationGrid->GetDataSize(); ++i)
    {
      pCorrelationGrid->SmearPoint(i);

    }

    return pScanMatcher;
  }

  /**
   * Match given scan against set of scans
   * @param pScan scan being scan-matched
   * @param rBaseScans set of scans whose points will mark cells in grid as being occupied
   * @param rMean output parameter of mean (best pose) of match
   * @param rCovariance output parameter of covariance of match
   * @param doPenalize whether to penalize matches further from the search center
   * @param doRefineMatch whether to do finer-grained matching if coarse match is good (default is true)
   * @return strength of response
   */
  double ScanMatcher::MatchScan(const LocalizedRangeScan& pScan, Pose2& rMean,
                                   Matrix3& rCovariance, bool doPenalize, bool doRefineMatch)
  {
    ////////////////////1///////////////////
    // set scan pose to be center of grid

    // 1. get scan position
    Pose2 scanPose = pScan.GetSensorPose();

    // scan has no readings; cannot do scan matching
    // best guess of pose is based off of adjusted odometer reading
    if (pScan.GetNumberOfRangeReadings() == 0)
    {
      rMean = scanPose;

      // maximum covariance
      rCovariance(0, 0) = MAX_VARIANCE;  // XX
      rCovariance(1, 1) = MAX_VARIANCE;  // YY
      rCovariance(2, 2) = 4 * math::Square(m_pMapper->m_pCoarseAngleResolution->GetValue());  // TH*TH

      return 0.0;
    }

    m_pCorrelationGrid->GetCoordinateConverter()->GetOffset();
    std::cout << m_pCorrelationGrid->GetCoordinateConverter()->GetOffset();

    // compute how far to search in each direction
    Vector2<double> searchDimensions(m_pSearchSpaceProbs->GetWidth(), m_pSearchSpaceProbs->GetHeight());

    Vector2<double> coarseSearchOffset(0.5 * (searchDimensions.GetX() - 1) * m_pCorrelationGrid->GetResolution(),
                                          0.5 * (searchDimensions.GetY() - 1) * m_pCorrelationGrid->GetResolution());

    // a coarse search only checks half the cells in each dimension
    Vector2<double> coarseSearchResolution( m_pCorrelationGrid->GetResolution(),
                                                m_pCorrelationGrid->GetResolution());

    // actual scan-matching
    double bestResponse = CorrelateScan(pScan, scanPose, coarseSearchOffset, coarseSearchResolution,
                                           m_pMapper->m_pCoarseSearchAngleOffset->GetValue(),
                                           m_pMapper->m_pCoarseAngleResolution->GetValue(),
                                           doPenalize, rMean, rCovariance, false);


// #ifdef KARTO_DEBUG
    // std::cout << "  BEST POSE = " << rMean << " BEST RESPONSE = " << bestResponse << ",  VARIANCE = "
              // << rCovariance(0, 0) << ", " << rCovariance(1, 1) << std::endl;
// #endif
    assert(math::InRange(rMean.GetHeading(), -PI, PI));

    return bestResponse;
  }

  /**
   * Finds the best pose for the scan centering the search in the correlation grid
   * at the given pose and search in the space by the vector and angular offsets
   * in increments of the given resolutions
   * @param rScan scan to match against correlation grid
   * @param rSearchCenter the center of the search space
   * @param rSearchSpaceOffset searches poses in the area offset by this vector around search center
   * @param rSearchSpaceResolution how fine a granularity to search in the search space
   * @param searchAngleOffset searches poses in the angles offset by this angle around search center
   * @param searchAngleResolution how fine a granularity to search in the angular search space
   * @param doPenalize whether to penalize matches further from the search center
   * @param rMean output parameter of mean (best pose) of match
   * @param rCovariance output parameter of covariance of match
   * @param doingFineMatch whether to do a finer search after coarse search
   * @return strength of response
   */
  double ScanMatcher::CorrelateScan(const LocalizedRangeScan& pScan, const Pose2& rSearchCenter,
                                       const Vector2<double>& rSearchSpaceOffset,
                                       const Vector2<double>& rSearchSpaceResolution,
                                       double searchAngleOffset, double searchAngleResolution,
                                       bool doPenalize, Pose2& rMean, Matrix3& rCovariance, bool doingFineMatch)
  {
    assert(searchAngleResolution != 0.0);

    // setup lookup arrays
    m_pGridLookup->ComputeOffsets(pScan, rSearchCenter.GetHeading(), searchAngleOffset, searchAngleResolution);

    // calculate position arrays

    std::vector<double> xPoses;
    unsigned int nX = static_cast<unsigned int>(math::Round(rSearchSpaceOffset.GetX() *
                                          2.0 / rSearchSpaceResolution.GetX()) + 1);

    double startX = -rSearchSpaceOffset.GetX();
    for (unsigned int xIndex = 0; xIndex < nX; xIndex++)
    {
      xPoses.push_back(startX + xIndex * rSearchSpaceResolution.GetX());
    }

    std::vector<double> yPoses;
    unsigned int nY = static_cast
    <unsigned int>(math::Round(rSearchSpaceOffset.GetY() *
                                          2.0 / rSearchSpaceResolution.GetY()) + 1);

    double startY = -rSearchSpaceOffset.GetY();
    for (unsigned int yIndex = 0; yIndex < nY; yIndex++)
    {
      yPoses.push_back(startY + yIndex * rSearchSpaceResolution.GetY());
    }
    // assert(math::DoubleEqual(yPoses.back(), -startY));

    // calculate pose response array size
    unsigned int nAngles = static_cast<unsigned int>(math::Round(searchAngleOffset * 2.0 / searchAngleResolution) + 1);

    unsigned int poseResponseSize = static_cast<unsigned int>(xPoses.size() * yPoses.size() * nAngles);

    // allocate array
    std::vector<std::pair<double, Pose2>> pPoseResponse;


    Vector2<int> startGridPoint = m_pCorrelationGrid->WorldToGrid(Vector2<double>(rSearchCenter.GetX()
                                                                        + startX, rSearchCenter.GetY() + startY),true);

    // std::cout<<"startGridPoint_World_X "<< rSearchCenter.GetX() + startX << std::endl;
    // std::cout<<"startGridPoint_World_Y "<< rSearchCenter.GetY() + startY << std::endl;

    // std::cout<<"startGridPoint_Grid_X "<< startGridPoint.GetX() << std::endl;
    // std::cout<<"startGridPoint_Grid_Y "<< startGridPoint.GetY() << std::endl;

    unsigned int poseResponseCounter = 0;
    forEachAs(std::vector<double>, &yPoses, yIter)
    {
      double y = *yIter;
      double newPositionY = rSearchCenter.GetY() + y;

      // std::cout<<"rSearchCenter_Y "<< rSearchCenter.GetY() << std::endl;

      double squareY = math::Square(y);

      forEachAs(std::vector<double>, &xPoses, xIter)
      {
        double x = *xIter;
        double newPositionX = rSearchCenter.GetX() + x;
        double squareX = math::Square(x);

        Vector2<int> gridPoint = m_pCorrelationGrid->WorldToGrid(Vector2<double>(newPositionX, newPositionY),true);
        int gridIndex = m_pCorrelationGrid->GridIndex(gridPoint);

        if(!math::IsUpTo(gridIndex, m_pCorrelationGrid->GetDataSize())){
          continue;
        }

        //assert(gridIndex >= 0);

        double angle = 0.0;
        double startAngle = rSearchCenter.GetHeading() - searchAngleOffset;
        for (unsigned int angleIndex = 0; angleIndex < nAngles; angleIndex++)
        {
          angle = startAngle + angleIndex * searchAngleResolution;

          double response = GetResponse(angleIndex, gridIndex);

          // store response and pose
          pPoseResponse.push_back(std::pair<double, Pose2>(response, Pose2(newPositionX, newPositionY,
                                                                            math::NormalizeAngle(angle))));
          poseResponseCounter++;
        }

        assert(poseResponseCounter == pPoseResponse.size());

        // assert(math::DoubleEqual(angle, rSearchCenter.GetHeading() + searchAngleOffset));
      }
    }

    assert(poseResponseSize == poseResponseCounter);

    // find value of best response (in [0; 1])
    double bestResponse = -1;
    int index_max_response = 0;
    for (unsigned int i = 0; i < poseResponseSize; i++)
    {
      // bestResponse = math::Maximum(bestResponse, pPoseResponse[i].first);

      if(bestResponse < pPoseResponse[i].first)
      {
        bestResponse = pPoseResponse[i].first;
        index_max_response  = i;
      }
      
    }

    //通过得分最大值位置之间的坐标差的欧式距离，利用正态分布概率密度函数求出权重计算机器人当前位置匹配得分
    double target_x = pPoseResponse[index_max_response].second.GetX();
    double target_y = pPoseResponse[index_max_response].second.GetY();

    double delta_x = target_x - rSearchCenter.GetX();
    double delta_y = target_y - rSearchCenter.GetY();

    double l2_norm = std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2));
    std::cout << "l2_norm = " << l2_norm << std::endl;

    double sigma = 0.2;
    //正态分布概率密度函数，N(0,sigma)
    double ratio = 1/(sqrt(2 * 3.14) * sigma) * exp(-std::pow(l2_norm, 2)  / (2 * sigma * sigma));
    double ratio_max =  1/(sqrt(2 * 3.14) * sigma);

    double ratio_target = ratio / ratio_max;

    std::cout << "bestResponse = " << bestResponse << std::endl;
    std::cout << "ratio_target = " << ratio_target << std::endl;

    bestResponse *= ratio_target;

    if (bestResponse > 1.0)
    {
      bestResponse = 1.0;
    }

    assert(math::InRange(bestResponse, 0.0, 1.0));
    // assert(math::InRange(rMean.GetHeading(), -PI, PI));

    return bestResponse;
  }

  

  /**
   * Get response at given position for given rotation (only look up valid points)
   * @param angleIndex
   * @param gridPositionIndex
   * @return response
   */
  double ScanMatcher::GetResponse(unsigned int angleIndex, int gridPositionIndex) const
  {
    double response = 0.0;
    // // // add up value for each point
    uint8_t* pByte = m_pCorrelationGrid->GetDataPointer() + gridPositionIndex;

    assert(pByte != NULL);

    const LookupArray* pOffsets = m_pGridLookup->GetLookupArray(angleIndex);

    // // get number of points in offset list
    unsigned int nPoints = pOffsets->GetSize();
    // std::cout << "pOffsets size = " << pOffsets->GetSize() << std::endl;
    if (nPoints == 0)
    {
      return response;
    }

    // calculate response
    std::vector<double> response_set;  
    std::map<int, double> response_map;
    const std::vector<unsigned int>& seg_counter = m_pGridLookup->GetSegCounter();

    // const std::vector<double>& slope_tmp = m_pGridLookup->GetSlope_set();

    unsigned int seed = 0;
    for(int i = 0; i < seg_counter.size(); ++i)
    {
      double tmp_response = 0;

      for(int j = 0; j < seg_counter[i]; ++j)
      {
        int index_offset = (*pOffsets)[seed];
        int pointGridIndex = gridPositionIndex + index_offset;
        if (!math::IsUpTo(pointGridIndex, m_pCorrelationGrid->GetDataSize()) || index_offset == INVALID_SCAN || pointGridIndex <= 0)
        {
          seed++;
          continue;
        }

        // uses index offsets to efficiently find location of point in the grid
        tmp_response += pByte[index_offset];
        seed++;
      }
      tmp_response /= (seg_counter[i] * GridStates_Occupied);
      response_set.push_back(tmp_response);
      response_map.insert(std::pair<int, double>(seg_counter[i], tmp_response));
    }

    // std::sort(response_set.begin(), response_set.end(), greater<double>());

    // uint8_t cal_size;
    // if(0 < response_set.size() < 5)
    // {
    //   cal_size = response_set.size();
    // }
    // else if(response_set.size() >= 5)
    // {
    //   cal_size = response_set.size() / 2;
    // }
    // else
    // {
    //   cal_size = 0;
    // }
    
    // double rr = 0;

    // if(cal_size > 0)
    // {
    //   for(int i = 0; i < cal_size; ++i)
    //   {
    //     rr += response_set[i];
    //   }
    //   response = rr / cal_size;
    // }
    
    //区块点数多的权重大，区块点数少的权重小
    std::map<int, double>::iterator it = response_map.begin();
    int num_sum =  0;
    for(; it != response_map.end(); ++it)
    {
      num_sum += it->first; 
    }

    it = response_map.begin();

    for(; it != response_map.end(); ++it)
    {
      response += it->second * it->first / num_sum; 
      std::cout << "seg size: " << it->first << '\t' << "response: " << it->second << std::endl;
    }

    std::cout<< "The score is = " << response << std::endl;
    
    // cal_size = response_map.size();
    
    // //判断得分最大的前几个之中是否有直线斜率接近相等的（过滤出来的直线有平行的），减小平行直线的权重，重点判断相交线段的得分
    // //相交直线存在情况下，平行直线的权重高，没有相交直线时，平行直线的权重低
    // std::vector<double> score_set;
    // for(int i = 0; i < cal_size; ++i, ++it)
    // {
    //   double tmp_previous;

    //   if(fabs(it->second - tmp_previous) < 0.3)
    //   {
    //     //mark当前直线，作为对比

    //   }else
    //   {        
    //   }
    //   //判断相交
    // }

    assert(fabs(response) <= 1.0);
    return response;
  }
}//namespace karto

