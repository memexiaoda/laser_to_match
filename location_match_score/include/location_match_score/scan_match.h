#ifndef CORRELATION_SCAN_MATCH_H
#define CORRELATION_SCAN_MATCH_H


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

#include "Parameter.h"
#include "Grid.h"
#include "data_structure.h"
#include "laser_scan.h"
#include "Matrix.h"

namespace karto
{
  #define MAX_VARIANCE            500.0
  #define DISTANCE_PENALTY_GAIN   0.2
  #define ANGLE_PENALTY_GAIN      0.2

/**
   * Scan matcher
   */
  class KARTO_EXPORT ScanMatcher
  {
  public:
    /**
     * Destructor
     */
    virtual ~ScanMatcher();

  public:
    /**
     * Create a scan matcher with the given parameters
     */
    static ScanMatcher* Create(Mapper* pMapper,
                               double searchSize,
                               double resolution,
                               double smearDeviation,
                               double rangeThreshold,
                               CorrelationGrid* pCorrelationGrid);

    double MatchScan(const LocalizedRangeScan& pScan,
                        Pose2& rMean, Matrix3& rCovariance,
                        bool doPenalize = true,
                        bool doRefineMatch = true);

    double CorrelateScan(const LocalizedRangeScan& pScan,
                            const Pose2& rSearchCenter,
                            const Vector2<double>& rSearchSpaceOffset,
                            const Vector2<double>& rSearchSpaceResolution,
                            double searchAngleOffset,
                            double searchAngleResolution,
                            bool doPenalize,
                            Pose2& rMean,
                            Matrix3& rCovariance,
                            bool doingFineMatch);


    // void ComputePositionalCovariance(const Pose2& rBestPose,
    //                                  double bestResponse,
    //                                  const Pose2& rSearchCenter,
    //                                  const Vector2<double>& rSearchSpaceOffset,
    //        location_match_score                          const Vector2<double>& rSearchSpaceResolution,
    //                                  double searchAngleResolution,
    //                                  Matrix3& rCovariance);

    // void ComputeAngularCovariance(const Pose2& rBestPose,
    //                               double bestResponse,
    //                               const Pose2& rSearchCenter,
    //                               double searchAngleOffset,
    //                               double searchAngleResolution,
    //                               Matrix3& rCovariance);

    /**
     * Gets the correlation grid data (for debugging)
     * @return correlation grid
     */
    inline CorrelationGrid* GetCorrelationGrid() const
    {
      return m_pCorrelationGrid;
    }

  private:
    /**
     * Marks cells where scans' points hit as being occupied
     * @param rScans scans whose points will mark cells in grid as being occupied
     * @param viewPoint do not add points that belong to scans "opposite" the view point
     */
    /**
     * Get response at given position for given rotation (only look up valid points)
     * @param angleIndex
     * @param gridPositionIndex
     * @return response
     */
    double GetResponse(unsigned int angleIndex, int gridPositionIndex) const;

  protected:
    /**
     * Default constructor
     */
    ScanMatcher(Mapper* pMapper)
      : m_pMapper(pMapper)
      , m_pCorrelationGrid(NULL)
      , m_pSearchSpaceProbs(NULL)
      , m_pGridLookup(NULL)
    {
    }

  private:
    Mapper* m_pMapper;
    CorrelationGrid* m_pCorrelationGrid;
    Grid<double>* m_pSearchSpaceProbs;
    GridIndexLookup<uint8_t>* m_pGridLookup;
  };  // ScanMatcher
}//namespace karto

#endif