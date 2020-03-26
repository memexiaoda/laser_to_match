/*
 * Copyright 2010 SRI International
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPEN_KARTO_MATH_H
#define OPEN_KARTO_MATH_H

#include <assert.h>
#include <math.h>
#include <limits>

#include "Types.h"

namespace karto
{
  /**
   * Platform independent pi definitions
   */
  const double PI         =  3.14159265358979323846;  // The value of PI
  const double Two_PI        =  6.28318530717958647692;  // 2 * PI
  const double PI_2       =  1.57079632679489661923;  // PI / 2
  const double PI_180     =  0.01745329251994329577;  // PI / 180
  const double Rad_2_Deg     = 57.29577951308232087685;  // 180 / PI

  /**
   * Lets define a small number!
   */
  const double TOLERANCE  = 1e-05;

  /**
     * Lets define max value of int (int32_t) to use it to mark invalid scans
     */

  const int INVALID_SCAN = std::numeric_limits<int>::max();

  namespace math
  {
    /**
     * Converts degrees into radians
     * @param degrees
     * @return radian equivalent of degrees
     */
    inline double DegreesToRadians(double degrees)
    {
      return degrees * PI_180;
    }

    /**
     * Converts radians into degrees
     * @param radians
     * @return degree equivalent of radians
     */
    inline double RadiansToDegrees(double radians)
    {
      return radians * Rad_2_Deg;
    }

    /**
     * Square function
     * @param value
     * @return square of value
     */
    template<typename T>
    inline T Square(T value)
    {
      return (value * value);
    }

    /**
     * Round function
     * @param value
     * @return rounds value to the nearest whole number (as double)
     */
    inline double Round(double value)
    {
      return value >= 0.0 ? floor(value + 0.5) : ceil(value - 0.5);
    }

    /**
     * Binary minimum function
     * @param value1
     * @param value2
     * @return the lesser of value1 and value2
     */
    template<typename T>
    inline const T& Minimum(const T& value1, const T& value2)
    {
      return value1 < value2 ? value1 : value2;
    }

    /**
     * Binary maximum function
     * @param value1
     * @param value2
     * @return the greater of value1 and value2
     */
    template<typename T>
    inline const T& Maximum(const T& value1, const T& value2)
    {
      return value1 > value2 ? value1 : value2;
    }

    /**
     * Clips a number to the specified minimum and maximum values.
     * @param n number to be clipped
     * @param minValue minimum value
     * @param maxValue maximum value
     * @return the clipped value
     */
    template<typename T>
    inline const T& Clip(const T& n, const T& minValue, const T& maxValue)
    {
      return Minimum(Maximum(n, minValue), maxValue);
    }

    /**
     * Checks whether two numbers are equal within a certain tolerance.
     * @param a
     * @param b
     * @return true if a and b differ by at most a certain tolerance.
     */
    inline bool DoubleEqual(double a, double b)
    {
      double delta = a - b;
      return delta < 0.0 ? delta >= -TOLERANCE : delta <= TOLERANCE;
    }

    /**
     * Checks whether value is in the range [0;maximum)
     * @param value
     * @param maximum
     */
    template<typename T>
    inline bool IsUpTo(const T& value, const T& maximum)
    {
      return (value >= 0 && value < maximum);
    }

    /**
     * Checks whether value is in the range [0;maximum)
     * Specialized version for unsigned int (unsigned int)
     * @param value
     * @param maximum
     */
    template<>
    inline bool IsUpTo<unsigned int>(const unsigned int& value, const unsigned int& maximum)
    {
      return (value < maximum);
    }


    /**
     * Checks whether value is in the range [a;b]
     * @param value
     * @param a
     * @param b
     */
    template<typename T>
    inline bool InRange(const T& value, const T& a, const T& b)
    {
      return (value >= a && value <= b);
    }

    /**
     * Normalizes angle to be in the range of [-pi, pi]
     * @param angle to be normalized
     * @return normalized angle
     */
    inline double NormalizeAngle(double angle)
    {
      while (angle < -PI)
      {
        if (angle < -Two_PI)
        {
          angle += (unsigned int)(angle / -Two_PI) * Two_PI;
        }
        else
        {
          angle += Two_PI;
        }
      }

      while (angle > PI)
      {
        if (angle > Two_PI)
        {
          angle -= (unsigned int)(angle / Two_PI) * Two_PI;
        }
        else
        {
          angle -= Two_PI;
        }
      }

      assert(math::InRange(angle, -PI, PI));

      return angle;
    }

    /**
     * Returns an equivalent angle to the first parameter such that the difference
     * when the second parameter is subtracted from this new value is an angle
     * in the normalized range of [-pi, pi], i.e. abs(minuend - subtrahend) <= pi.
     * @param minuend
     * @param subtrahend
     * @return normalized angle
     */
    inline double NormalizeAngleDifference(double minuend, double subtrahend)
    {
      while (minuend - subtrahend < -PI)
      {
        minuend += Two_PI;
      }

      while (minuend - subtrahend > PI)
      {
        minuend -= Two_PI;
      }

      return minuend;
    }

    /**
     * Align a value to the alignValue.
     * The alignValue should be the power of two (2, 4, 8, 16, 32 and so on)
     * @param value
     * @param alignValue
     * @return aligned value
     */
    template<class T>
    inline T AlignValue(size_t value, size_t alignValue = 8)
    {
      return static_cast<T> ((value + (alignValue - 1)) & ~(alignValue - 1));
    }
  }  // namespace math

}  // namespace karto

#endif  // OPEN_KARTO_MATH_H
