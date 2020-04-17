#ifndef _LASER_RANGE_FINDER_H
#define _LASER_RANGE_FINDER_H

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
#include <memory>

#include "Parameter.h"
#include "Grid.h"
#include "data_structure.h"
#include "laser_scan.h"
#include "Matrix.h"


namespace karto
{
/**
   * The LaserRangeFinder defines a laser sensor that provides the pose offset position of a localized range scan relative to the robot.
   * The user can set an offset pose for the sensor relative to the robot coordinate system. If no value is provided
   * by the user, the sensor is set to be at the origin of the robot coordinate system.
   * The LaserRangeFinder contains parameters for physical laser sensor used by the mapper for scan matching
   * Also contains information about the maximum range of the sensor and provides a threshold
   * for limiting the range of readings.
   * The optimal value for the range threshold depends on the angular resolution of the scan and
   * the desired map resolution.  RangeThreshold should be set as large as possible while still
   * providing "solid" coverage between consecutive range readings.  The diagram below illustrates
   * the relationship between map resolution and the range threshold.
   */
  class KARTO_EXPORT LaserRangeFinder
  {

  public:
    /**
     * Destructor
     */
    ~LaserRangeFinder() {}

    /**
     * Gets this range finder sensor's offset
     * @return offset pose
     */
    inline const Pose2 GetOffsetPose() const
    {
      return m_pOffsetPose->GetValue();
    }


  /**
   * Sets this range finder sensor's offset
   * @param rPose
  */
    inline void SetOffsetPose(const Pose2 rPose)
    {
      m_pOffsetPose->SetValue(rPose);
    }

  public:
    /**
     * Gets this range finder sensor's minimum range
     * @return minimum range
     */
    inline double GetMinimumRange() const
    {
      return m_pMinimumRange->GetValue();
    }

    /**
     * Sets this range finder sensor's minimum range
     * @param minimumRange
     */
    inline void SetMinimumRange(double minimumRange)
    {
      m_pMinimumRange->SetValue(minimumRange);

      SetRangeThreshold(GetRangeThreshold());
    }

    /**
     * Gets this range finder sensor's maximum range
     * @return maximum range
     */
    inline double GetMaximumRange() const
    {
      return m_pMaximumRange->GetValue();
    }

    /**
     * Sets this range finder sensor's maximum range
     * @param maximumRange
     */
    inline void SetMaximumRange(double maximumRange)
    {
      m_pMaximumRange->SetValue(maximumRange);

      SetRangeThreshold(GetRangeThreshold());
    }

    /**
     * Gets the range threshold
     * @return range threshold
     */
    inline double GetRangeThreshold() const
    {
      return m_pRangeThreshold->GetValue();
    }

    /**
     * Sets the range threshold
     * @param rangeThreshold
     */
    inline void SetRangeThreshold(double rangeThreshold)
    {
      // make sure rangeThreshold is within laser range finder range
      m_pRangeThreshold->SetValue(math::Clip(rangeThreshold, GetMinimumRange(), GetMaximumRange()));

      if (math::DoubleEqual(GetRangeThreshold(), rangeThreshold) == false)
      {
       // std::cout << "Info: clipped range threshold to be within minimum and maximum range!" << std::endl;
      }
    }

    /**
     * Gets this range finder sensor's minimum angle
     * @return minimum angle
     */
    inline double GetMinimumAngle() const
    {
      return m_pMinimumAngle->GetValue();
    }

    /**
     * Sets this range finder sensor's minimum angle
     * @param minimumAngle
     */
    inline void SetMinimumAngle(double minimumAngle)
    {
      m_pMinimumAngle->SetValue(minimumAngle);

      Update();
    }

    /**
     * Gets this range finder sensor's maximum angle
     * @return maximum angle
     */
    inline double GetMaximumAngle() const
    {
      return m_pMaximumAngle->GetValue();
    }

    /**
     * Sets this range finder sensor's maximum angle
     * @param maximumAngle
     */
    inline void SetMaximumAngle(double maximumAngle)
    {
      m_pMaximumAngle->SetValue(maximumAngle);

      Update();
    }

    /**
     * Gets this range finder sensor's angular resolution
     * @return angular resolution
     */
    inline double GetAngularResolution() const
    {
      return m_pAngularResolution->GetValue();
    }

    /**
     * Sets this range finder sensor's angular resolution
     * @param angularResolution
     */
    inline void SetAngularResolution(double angularResolution)
    {
      
      
      m_pAngularResolution->SetValue(angularResolution);
      

      Update();
    }

    inline unsigned int GetNumberOfRangeReadings()
    {
      return m_NumberOfRangeReadings;

    }

    /**
     * Gets the number of range readings each localized range scan must contain to be a valid scan.
     * @return number of range readings
     */
    inline unsigned int GetNumberOfRangeReadings() const
    {
      return m_NumberOfRangeReadings;
    }

    virtual bool Validate()
    {
      Update();

      if (math::InRange(GetRangeThreshold(), GetMinimumRange(), GetMaximumRange()) == false)
      {
       // std::cout << "Please set range threshold to a value between ["
        //          << GetMinimumRange() << ";" << GetMaximumRange() << "]" << std::endl;
        return false;
      }

      return true;
    }



//    virtual bool Validate(SensorData* pSensorData);

  public:
    /**
     * Create a laser range finder of the given type and ID
     * @param type
     * @param rName name of sensor - if no name is specified default name will be assigned
     * @return laser range finder
     */
    static std::shared_ptr<LaserRangeFinder> CreateLaserRangeFinder()
    {
      std::shared_ptr<LaserRangeFinder> p_Laser(new LaserRangeFinder());
      return p_Laser;
    }

  private:
    /**
     * Constructs a LaserRangeFinder object with given ID
     */
    LaserRangeFinder()
      : m_NumberOfRangeReadings(0),
        m_pOffsetPose(new Parameter<Pose2>(Pose2())),
        m_pMinimumRange(new Parameter<double>(0.0)),
        m_pMaximumRange(new Parameter<double>(80.0)),
        m_pMinimumAngle(new Parameter<double>(math::DegreesToRadians(-90))),
        m_pMaximumAngle(new Parameter<double>(math::DegreesToRadians(90))),
        m_pAngularResolution(new Parameter<double>(math::DegreesToRadians(1.0))),
        m_pRangeThreshold(new Parameter<double>(30.0))
    {

    }

    /**
     * Set the number of range readings based on the minimum and
     * maximum angles of the sensor and the angular resolution
     */
    void Update()
    {
      m_NumberOfRangeReadings = static_cast<unsigned int>(math::Round((GetMaximumAngle() -
                                                                    GetMinimumAngle())
                                                                    / GetAngularResolution()) + 1);
    }

  private:
    LaserRangeFinder(const LaserRangeFinder&);
    const LaserRangeFinder& operator=(const LaserRangeFinder&);

  private:
    // sensor m_Parameters
    std::unique_ptr<Parameter<Pose2>> m_pOffsetPose;
    std::unique_ptr<Parameter<double>> m_pMinimumAngle;
    std::unique_ptr<Parameter<double>> m_pMaximumAngle;

    std::unique_ptr<Parameter<double>> m_pAngularResolution;

    std::unique_ptr<Parameter<double>> m_pMinimumRange;
    std::unique_ptr<Parameter<double>> m_pMaximumRange;

    std::unique_ptr<Parameter<double>> m_pRangeThreshold;

    unsigned int m_NumberOfRangeReadings;

    // static std::string LaserRangeFinderTypeNames[6];
  };  // LaserRangeFinder


}

  #endif