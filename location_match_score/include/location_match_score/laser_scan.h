#ifndef _LASER_SCAN_H
#define _LASER_SCAN_H

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
#include <boost/thread.hpp>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "data_structure.h"
#include "Pose.h"
#include "Transform.h"
#include "laser_range_finder.h"

#include <eigen3/Eigen/Dense>


namespace karto
{
  /**
   * Type declaration of range readings vector
   */
  typedef std::vector<Eigen::Vector2f> RangeReadingsVector;

 /**
   * LaserRangeScan representing the range readings from a laser range finder sensor.
   */
  class LaserRangeScan
  {

  public:
    /**
     * Constructs a scan from the given sensor with the given readings
     * @param rSensorName
     */
    LaserRangeScan()
      :m_NumberOfRangeReadings(0)
    {
    }

    /**
     * Constructs a scan from the given sensor with the given readings
     * @param rSensorName
     * @param rRangeReadings
     */
    LaserRangeScan(const RangeReadingsVector& rRangeReadings)
      :m_NumberOfRangeReadings(0)
    {

      SetRangeReadings(rRangeReadings);
    }

    /**
     * Destructor
     */
    virtual ~LaserRangeScan()
    {
      // std::cout<< "m_pRangeReadings release." <<std::endl;
    }

  public:
    /**
     * Gets the range readings of this scan
     * @return range readings of this scan
     */
    inline const std::vector<Eigen::Vector2f>* GetRangeReadings() const
    {
      return &m_RangeReadings;
    }

    inline RangeReadingsVector GetRangeReadingsVector() const
    {
      return m_RangeReadings;
    }

    /**
     * Sets the range readings for this scan
     * @param rRangeReadings
     */
    inline void SetRangeReadings(const RangeReadingsVector& rRangeReadings)
    {

      if (!rRangeReadings.empty())
      {
        
          // delete old readings
          m_RangeReadings.clear();

          // store size of array!
          m_NumberOfRangeReadings = static_cast<unsigned int>(rRangeReadings.size());
        
        // copy readings
        unsigned int index = 0;
        const_forEach(RangeReadingsVector, &rRangeReadings)
        {
          m_RangeReadings.push_back(*iter);
        }

        // std::cout << "copy readings number = "<< m_RangeReadings.size() << std::endl;
      }
      else
      {
        m_RangeReadings.clear();
      }
    }

    /**
     * Gets the number of range readings
     * @return number of range readings
     */
    inline unsigned int GetNumberOfRangeReadings() const
    {
      return m_NumberOfRangeReadings;
    }

  private:
    LaserRangeScan(const LaserRangeScan&);
    const LaserRangeScan& operator=(const LaserRangeScan&);

  private:
    std::vector<Eigen::Vector2f> m_RangeReadings;
    unsigned int m_NumberOfRangeReadings;
  };  // LaserRangeScan
  
  class LocalizedRangeScan : public LaserRangeScan
  {

  public:
    /**
     * Constructs a range scan from the given range finder with the given readings
     */
    LocalizedRangeScan()
      : m_IsDirty(true)
    {
    }

    /**
     * Destructor
     */
    virtual ~LocalizedRangeScan()
    {
    }

  private:
    mutable boost::shared_mutex m_Lock;

  public:
    /**
     * Gets the odometric pose of this scan
     * @return odometric pose of this scan
     */
    inline const Pose2& GetOdometricPose() const
    {
      return m_OdometricPose;
    }

    /**
     * Sets the odometric pose of this scan
     * @param rPose
     */
    inline void SetOdometricPose(const Pose2& rPose)
    {
      m_OdometricPose = rPose;
    }

    /**
     * Gets the (possibly corrected) robot pose at which this scan was taken.  The corrected robot pose of the scan
     * is usually set by an external module such as a localization or mapping module when it is determined
     * that the original pose was incorrect.  The external module will set the correct pose based on
     * additional sensor data and any context information it has.  If the pose has not been corrected,
     * a call to this method returns the same pose as GetOdometricPose().
     * @return corrected pose
     */
    inline const Pose2& GetCorrectedPose() const
    {
      return m_CorrectedPose;
    }

    /**
     * Moves the scan by moving the robot pose to the given location.
     * @param rPose new pose of the robot of this scan
     */
    inline void SetCorrectedPose(const Pose2& rPose)
    {
      m_CorrectedPose = rPose;

      m_IsDirty = true;
    }

    inline void SetLaserRangeFinder(std::shared_ptr<LaserRangeFinder> pLaser)
    {
       m_pLaserRangeFinder = pLaser;
    }

    /**
     * Computes the position of the sensor
     * @return scan pose
     */
     inline Pose2 GetSensorPose() const
     {
       //return GetSensorAt(m_CorrectedPose);
        return m_CorrectedPose;
     }

    /**
     * Get point readings in local coordinates
     */
     inline const PointVectorDouble& GetPointReadings(bool wantFiltered = false) const
     {
        boost::shared_lock<boost::shared_mutex> lock(m_Lock);
        if (m_IsDirty)
        {
          // throw away constness and do an update!
          lock.unlock();
          boost::unique_lock<boost::shared_mutex> uniqueLock(m_Lock);
          const_cast<LocalizedRangeScan*>(this)->Update();
        }

        return m_PointReadings;
     }

     const PointVectorDouble& GetLocalPointReadings() const
     {
       return m_LocalPointReadings;
     }

  private:
    /**
     * Compute point readings based on range readings
     * Only range readings within [minimum range; range threshold] are returned
     */

    virtual void Update()
    {
      Pose2 scanPose = GetSensorPose();
      double heading = scanPose.GetHeading();
      // std::cout << "scanPose = " << scanPose << std::endl;

      m_PointReadings.clear();
      m_LocalPointReadings.clear();

    // compute point readings
       unsigned int beamNum = 0;
       const std::vector<Eigen::Vector2f>* pRangeReadings = GetRangeReadings();

      //  std::cout<<"number of RangeReadings = "<<m_pLaserRangeFinder->GetNumberOfRangeReadings()<<std::endl;

       for (unsigned int i = 0; i < GetRangeReadingsVector().size(); i++)
       {

         Eigen::Vector2f rangeReading = (*pRangeReadings)[i];

         Vector2<double> local_point;
         local_point.SetX(rangeReading.x());
         local_point.SetY(rangeReading.y());
         m_LocalPointReadings.push_back(local_point);

         Vector2<double> point;

         point.SetX(scanPose.GetX() + rangeReading.x() * cos(heading) - rangeReading.y() * sin(heading));//全局坐标系下的位姿
         point.SetY(scanPose.GetY() + rangeReading.x() * sin(heading) - rangeReading.y() * cos(heading));//全局坐标系下的位姿

         m_PointReadings.push_back(point);
       }

      // std::cout<< "m_PointReadings size is: "<<m_PointReadings.size()<<std::endl;

      m_IsDirty = false;

    }

  private:
    LocalizedRangeScan(const LocalizedRangeScan&);
    const LocalizedRangeScan& operator=(const LocalizedRangeScan&);

  private:
    /**
     * Odometric pose of robot
     */
    Pose2 m_OdometricPose;

    /**
     * Corrected pose of robot calculated by mapper (or localizer)
     */
    Pose2 m_CorrectedPose;

  protected:

    std::shared_ptr<LaserRangeFinder> m_pLaserRangeFinder;
    /**
     * Average of all the point readings
     */
    Pose2 m_BarycenterPose;

    /**
     * Vector of point readings
     */
    PointVectorDouble m_PointReadings;

    PointVectorDouble m_LocalPointReadings;

    /**
     * Vector of unfiltered point readings
     */
    PointVectorDouble m_UnfilteredPointReadings;

    /**
     * Bounding box of localized range scan
     */
    BoundingBox2 m_BoundingBox;

    /**
     * Internal flag used to update point readings, barycenter and bounding box
     */
    bool m_IsDirty;
  };  // LocalizedRangeScan

   /**
   * Type declaration of LocalizedRangeScan vector
   */
  typedef std::vector<LocalizedRangeScan*> LocalizedRangeScanVector;

}//namespace karto

#endif
