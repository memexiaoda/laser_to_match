#ifndef _GRID_H
#define _GRID_H

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

#include "data_structure.h"
#include "laser_scan.h"
#include "exception.h"
#include "segment.h"
// #include "fitLine.h"
#include <yaml-cpp/yaml.h>

namespace karto
{

typedef enum
  {
    GridStates_Unknown = 0,
    GridStates_Occupied = 100,
    GridStates_Free = 0
  } GridStates;

/**
   * The CoordinateConverter class is used to convert coordinates between world and grid coordinates
   * In world coordinates 1.0 = 1 meter where 1 in grid coordinates = 1 pixel!
   * Default scale for coordinate converter is 20 that converters to 1 pixel = 0.05 meter
   */
  class CoordinateConverter
  {
  public:
    /**
     * Default constructor
     */
    CoordinateConverter()
      : m_Scale(20.0)
    {
    }

  public:
    /**
     * Scales the value
     * @param value
     * @return scaled value
     */
    inline double Transform(double value)
    {
      return value * m_Scale;
    }

    /**
     * Converts the point from world coordinates to grid coordinates
     * @param rWorld world coordinate
     * @param flipY
     * @return grid coordinate
     */
    inline Vector2<int> WorldToGrid(const Vector2<double>& rWorld, bool flipY = false) const
    { 
      double gridX = (rWorld.GetX() - m_Offset.GetX()) * m_Scale;
      double gridY = 0.0;

      if (flipY == false)
      {
        gridY = (rWorld.GetY() - m_Offset.GetY()) * m_Scale;
      }
      else
      {
        gridY = (m_Size.GetHeight() / m_Scale - rWorld.GetY() + m_Offset.GetY()) * m_Scale;
      }

      return Vector2<int>(static_cast<int>(math::Round(gridX)), static_cast<int>(math::Round(gridY)));
    }

    /**
     * Converts the point from grid coordinates to world coordinates
     * @param rGrid world coordinate
     * @param flipY
     * @return world coordinate
     */
    inline Vector2<double> GridToWorld(const Vector2<int>& rGrid, bool flipY = false) const
    {
      double worldX = m_Offset.GetX() + rGrid.GetX() / m_Scale;
      double worldY = 0.0;

      if (flipY == false)
      {
        worldY = m_Offset.GetY() + rGrid.GetY() / m_Scale;
      }
      else
      {
        worldY = m_Offset.GetY() + (m_Size.GetHeight() - rGrid.GetY()) / m_Scale;
      }

      return Vector2<double>(worldX, worldY);
    }

    /**
     * Gets the scale
     * @return scale
     */
    inline double GetScale() const
    {
      return m_Scale;
    }

    /**
     * Sets the scale
     * @param scale
     */
    inline void SetScale(double scale)
    {
      m_Scale = scale;
    }

    /**
     * Gets the offset
     * @return offset
     */
    inline const Vector2<double>& GetOffset() const
    {
      return m_Offset;
    }

    /**
     * Sets the offset
     * @param rOffset
     */
    inline void SetOffset(const Vector2<double>& rOffset)
    {
      m_Offset = rOffset;
    }

    /**
     * Sets the size
     * @param rSize
     */
    inline void SetSize(const Size2<int>& rSize)
    {
      m_Size = rSize;
    }

    /**
     * Gets the size
     * @return size
     */
    inline const Size2<int>& GetSize() const
    {
      return m_Size;
    }

    /**
     * Gets the resolution
     * @return resolution
     */
    inline double GetResolution() const
    {
      return 1.0 / m_Scale;
    }

    /**
     * Sets the resolution
     * @param resolution
     */
    inline void SetResolution(double resolution)
    {
      m_Scale = 1.0 / resolution;
    }

    /**
     * Gets the bounding box
     * @return bounding box
     */
    inline BoundingBox2 GetBoundingBox() const
    {
      BoundingBox2 box;

      double minX = GetOffset().GetX();
      double minY = GetOffset().GetY();
      double maxX = minX + GetSize().GetWidth() * GetResolution();
      double maxY = minY + GetSize().GetHeight() * GetResolution();

      box.SetMinimum(GetOffset());
      box.SetMaximum(Vector2<double>(maxX, maxY));
      return box;
    }

  private:
    Size2<int> m_Size;
    double m_Scale;

    Vector2<double> m_Offset;
  };  // CoordinateConverter



/**
   * Defines a grid class
   */
  template<typename T>
  class Grid
  {
  public:
    /**
     * Creates a grid of given size and resolution
     * @param width
     * @param height
     * @param resolution
     * @return grid pointer
     */
    static Grid* CreateGrid(int width, int height, double resolution)
    {
      Grid* pGrid = new Grid(width, height);

      pGrid->GetCoordinateConverter()->SetScale(1.0 / resolution);

      return pGrid;
    }

    /**
     * Destructor
     */
    virtual ~Grid()
    {
      delete [] m_pData;
      delete m_pCoordinateConverter;
      std::cout<< "Deconstruction Grid" <<std::endl;
    }

  public:
    /**
     * Clear out the grid data
     */
    void Clear()
    {
      memset(m_pData, 0, GetDataSize() * sizeof(T));
    }

    /**
     * Returns a clone of this grid
     * @return grid clone
     */
    Grid* Clone()
    {
      Grid* pGrid = CreateGrid(GetWidth(), GetHeight(), GetResolution());
      pGrid->GetCoordinateConverter()->SetOffset(GetCoordinateConverter()->GetOffset());

      memcpy(pGrid->GetDataPointer(), GetDataPointer(), GetDataSize());

      return pGrid;
    }

    /**
     * Resizes the grid (deletes all old data)
     * @param width
     * @param height
     */
    virtual void Resize(int width, int height)
    {
      m_Width = width;
      m_Height = height;
      // m_WidthStep = math::AlignValue<int>(width, 8);

      if (m_pData != NULL)
      {
        delete[] m_pData;
        m_pData = NULL;
      }

      try
      {
        m_pData = new T[GetDataSize()];

        if (m_pCoordinateConverter == NULL)
        {
          m_pCoordinateConverter = new CoordinateConverter();
        }

        m_pCoordinateConverter->SetSize(Size2<int>(width, height));
      }
      catch(...)
      {
        m_pData = NULL;

        m_Width = 0;
        m_Height = 0;
        // m_WidthStep = 0;
      }

      //Clear();
    }

    /**
     * Checks whether the given coordinates are valid grid indices
     * @param rGrid
     */
    inline bool IsValidGridIndex(const Vector2<int>& rGrid) const
    {
      return (math::IsUpTo(rGrid.GetX(), m_Width) && math::IsUpTo(rGrid.GetY(), m_Height));
    }

    /**
     * Gets the index into the data pointer of the given grid coordinate
     * @param rGrid
     * @param boundaryCheck default value is true
     * @return grid index
     */
    virtual int GridIndex(const Vector2<int>& rGrid, bool boundaryCheck = false) const
    {
      // std::cout << "rGrid X = " << rGrid.GetX() <<std::endl;
      // std::cout << "rGrid Y = " << rGrid.GetY() <<std::endl;

      if (boundaryCheck == true)
      {
        if (IsValidGridIndex(rGrid) == false)
        {
          std::cout << "calculate index!"<<std::endl;
          std::stringstream error;
          error << "Index " << rGrid << " out of range.  Index must be between [0; "
                << m_Width << ") and [0; " << m_Height << ")";
          throw Exception(error.str());
        }
      }

      // int index = rGrid.GetX() + (rGrid.GetY() * m_WidthStep);

      int index = rGrid.GetX() + (rGrid.GetY() * m_Width);

      if (boundaryCheck == true)
      {
        assert(math::IsUpTo(index, GetDataSize()));
      }

      return index;
    }

    /**
     * Gets the grid coordinate from an index
     * @param index
     * @return grid coordinate
     */
    Vector2<int> IndexToGrid(int index) const
    {
      Vector2<int> grid;

      // grid.SetY(index / m_WidthStep);
      //grid.SetX(index - grid.GetY() * m_WidthStep);
      grid.SetY(index / m_Width);
      grid.SetX(index - grid.GetY() * m_Width);

      return grid;
    }

    /**
     * Converts the point from world coordinates to grid coordinates
     * @param rWorld world coordinate
     * @param flipY
     * @return grid coordinate
     */
    inline Vector2<int> WorldToGrid(const Vector2<double>& rWorld, bool flipY = false) const
    {
      return GetCoordinateConverter()->WorldToGrid(rWorld, flipY);
    }

    /**
     * Converts the point from grid coordinates to world coordinates
     * @param rGrid world coordinate
     * @param flipY
     * @return world coordinate
     */
    inline Vector2<double> GridToWorld(const Vector2<int>& rGrid, bool flipY = false) const
    {
      return GetCoordinateConverter()->GridToWorld(rGrid, flipY);
    }

    /**
     * Gets pointer to data at given grid coordinate
     * @param rGrid grid coordinate
     * @return grid point
     */
    T* GetDataPointer(const Vector2<int>& rGrid)
    {
      int index = GridIndex(rGrid, false);
      return m_pData + index;
    }

    /**
     * Gets pointer to data at given grid coordinate
     * @param rGrid grid coordinate
     * @return grid point
     */
    T* GetDataPointer(const Vector2<int>& rGrid) const
    {
      int index = GridIndex(rGrid, false);
      return m_pData + index;
    }

    /**
     * Gets the width of the grid
     * @return width of the grid
     */
    inline int GetWidth() const
    {
      return m_Width;
    };

    /**
     * Gets the height of the grid
     * @return height of the grid
     */
    inline int GetHeight() const
    {
      return m_Height;
    };

    /**
     * Get the size as a Size2<int>
     * @return size of the grid
     */
    inline const Size2<int> GetSize() const
    {
      return Size2<int>(m_Width, m_Height);
    }

    /**
     * Gets the grid data pointer
     * @return data pointer
     */
    inline T* GetDataPointer()
    {
      return m_pData;
    }

    /**
     * Gets const grid data pointer
     * @return data pointer
     */
    inline T* GetDataPointer() const
    {
      return m_pData;
    }

    /**
     * Gets the allocated grid size in bytes
     * @return data size
     */
    inline int GetDataSize() const
    {
      // return m_WidthStep * m_Height;

      return m_Width * m_Height;
    }

    /**
     * Get value at given grid coordinate
     * @param rGrid grid coordinate
     * @return value
     */
    inline T GetValue(const Vector2<int>& rGrid) const
    {
      int index = GridIndex(rGrid);
      return m_pData[index];
    }

    /**
     * Gets the coordinate converter for this grid
     * @return coordinate converter
     */
    inline CoordinateConverter* GetCoordinateConverter() const
    {
      return m_pCoordinateConverter;
    }

    /**
     * Gets the resolution
     * @return resolution
     */
    inline double GetResolution() const
    {
      return GetCoordinateConverter()->GetResolution();
    }

    /**
     * Gets the grids bounding box
     * @return bounding box
     */
    inline BoundingBox2 GetBoundingBox() const
    {
      return GetCoordinateConverter()->GetBoundingBox();
    }

  protected:
    /**
     * Constructs grid of given size
     * @param width
     * @param height
     */
    Grid(int width, int height)
      : m_pData(NULL)
      , m_pCoordinateConverter(NULL)
    {
      Resize(width, height);
    }

  private:
    int m_Width;       // width of grid
    int m_Height;      // height of grid
    T* m_pData;        // grid data

    // coordinate converter to convert between world coordinates and grid coordinates
    CoordinateConverter* m_pCoordinateConverter;
  };  // Grid



  class CorrelationGrid;
  
  /**
   * Create lookup tables for point readings at varying angles in grid.
   * For each angle, grid indexes are calculated for each range reading.
   * This is to speed up finding best angle/position for a localized range scan
   *
   * Used heavily in mapper and localizer.
   *
   * In the localizer, this is a huge speed up for calculating possible position.  For each particle,
   * a probability is calculated.  The range scan is the same, but all grid indexes at all possible angles are
   * calculated.  So when calculating the particle probability at a specific angle, the index table is used
   * to look up probability in probability grid!
   *
   */
  
  /**
   * An array that can be resized as long as the size
   * does not exceed the initial capacity
   */
  class LookupArray
  {
  public:
    /**
     * Constructs lookup array
     */
    LookupArray()
      : m_pArray(NULL)
      , m_Capacity(0)
      , m_Size(0)
    {
    }

    /**
     * Destructor
     */
    virtual ~LookupArray()
    {
      assert(m_pArray != NULL);

      delete[] m_pArray;
      m_pArray = NULL;
    }

  public:
    /**
     * Clear array
     */
    void Clear()
    {
      memset(m_pArray, 0, sizeof(int) * m_Capacity);
    }

    /**
     * Gets size of array
     * @return array size
     */
    unsigned int GetSize() const
    {
      return m_Size;
    }

    /**
     * Sets size of array (resize if not big enough)
     * @param size
     */
    void SetSize(unsigned int size)
    {
      assert(size != 0);

      if (size > m_Capacity)
      {
        if (m_pArray != NULL)
        {
          delete [] m_pArray;
        }
        m_Capacity = size;
        m_pArray = new int[m_Capacity];
      }

      m_Size = size;
    }

    /**
     * Gets reference to value at given index
     * @param index
     * @return reference to value at index
     */
    inline int& operator [] (unsigned int index)
    {
      assert(index < m_Size);

      return m_pArray[index];
    }

    /**
     * Gets value at given index
     * @param index
     * @return value at index
     */
    inline int operator [] (unsigned int index) const
    {
      // std::cout << "index = " << index << std::endl;
      // std::cout << "m_Size = " << m_Size << std::endl;
      assert(index < m_Size);

      return m_pArray[index];
    }

    /**
     * Gets array pointer
     * @return array pointer
     */
    inline int* GetArrayPointer()
    {
      return m_pArray;
    }

    /**
     * Gets array pointer
     * @return array pointer
     */
    inline int* GetArrayPointer() const
    {
      return m_pArray;
    }

  private:
    int* m_pArray;
    unsigned int m_Capacity;
    unsigned int m_Size;
  };  // LookupArray

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Create lookup tables for point readings at varying angles in grid.
   * For each angle, grid indexes are calculated for each range reading.
   * This is to speed up finding best angle/position for a localized range scan
   *
   * Used heavily in mapper and localizer.
   *
   * In the localizer, this is a huge speed up for calculating possible position.  For each particle,
   * a probability is calculated.  The range scan is the same, but all grid indexes at all possible angles are
   * calculated.  So when calculating the particle probability at a specific angle, the index table is used
   * to look up probability in probability grid!
   *
   */
  template<typename T>
  class GridIndexLookup
  {
  public:
    /**
     * Construct a GridIndexLookup with a grid
     * @param pGrid
     */
    GridIndexLookup(Grid<T>* pGrid)
      : m_pGrid(pGrid)
      , m_Capacity(0)
      , m_Size(0)
      , m_ppLookupArray(NULL)
      , m_dis_threshold(0.3)
      , m_batch_size(30)
    {
      LoadSegParameters();
    }

    /**
     * Destructor
     */
    virtual ~GridIndexLookup()
    {
      DestroyArrays();
    }

  public:
    /**
     * Gets the lookup array for a particular angle index
     * @param index
     * @return lookup array
     */
    const LookupArray* GetLookupArray(unsigned int index) const
    {

      assert(math::IsUpTo(index, m_Size));

      return m_ppLookupArray[index];
    }

    /**
     * Get angles
     * @return std::vector<double>& angles
     */
    const std::vector<double>& GetAngles() const
    {
      return m_Angles;
    }

    const std::vector<unsigned int>& GetSegCounter() const
    {
      return seg_counter;
    }

    const std::vector<double>& GetSlope_set() const
    {
      return slope_set;
    }

    const std::vector<double>& GetDisSet() const
    {
      return dis_set;
    }

    void LoadSegParameters()
    {
      YAML::Node config = YAML::LoadFile("/home/lei/scan_to_match/src/location_match_score/config/segment_params.yaml");
	    m_dis_threshold = config["dis_threshold"].as<float>();
	    m_batch_size = config["batch_size"].as<int>();
    }

    /**
     * Compute lookup table of the points of the given scan for the given angular space
     * @param pScan the scan
     * @param angleCenter
     * @param angleOffset computes lookup arrays for the angles within this offset around angleStart
     * @param angleResolution how fine a granularity to compute lookup arrays in the angular space
     */
    void ComputeOffsets(const LocalizedRangeScan& pScan,
                        double angleCenter,
                        double angleOffset,
                        double angleResolution)
    {
      assert(angleOffset != 0.0);
      assert(angleResolution != 0.0);

      unsigned int nAngles = static_cast<unsigned int>(math::Round(angleOffset * 2.0 / angleResolution) + 1);
      std::cout << "nAngles = " << nAngles << std::endl;
      SetSize(nAngles);

      // convert points into local coordinates of scan pose

      const PointVectorDouble& rPointReadings = pScan.GetPointReadings();

      // compute transform to scan pose
      Transform transform(pScan.GetSensorPose());      

      Pose2Vector localPoints;
      const PointVectorDouble& rLocalPointReadings = pScan.GetLocalPointReadings();
      const_forEach(PointVectorDouble, &rLocalPointReadings)
      {
        Pose2 vec = Pose2(*iter, 0.0);
        localPoints.push_back(vec);
      }
      std::cout << "localPoints size = " << localPoints.size() << std::endl;
      
      //过滤数据，将数据通过聚类方法分隔为小块
      seg_counter.clear();
      std::vector<Pose2Vector> segment_set;
      segment_set = segment(localPoints, m_dis_threshold, m_batch_size);

      Pose2Vector localPoints_filter;

      //这里增加一个判断是避免分割函数没有数据时，会使得程序异常终止
      if(segment_set.size() == 0)
      {
        segment_set.push_back(localPoints);
      }

      for(int i = 0; i < segment_set.size(); ++i)
      {
        unsigned int ss = segment_set[i].size();

        for(int j = 0; j < ss; ++j)
        {
          localPoints_filter.push_back((segment_set[i])[j]);
        }
        
        seg_counter.push_back(ss);

        std::cout << "seg " << i << "size = " << ss << std::endl;
      }

      std::cout << "localPoints_filter size() = " << localPoints_filter.size() << std::endl;

      const Vector2<double>& rGridOffset = m_pGrid->GetCoordinateConverter()->GetOffset();
       double cosine = cos(angleCenter);
      double sine = sin(angleCenter);

      slope_set.clear();

      dis_set.clear();

      //记录激光点距离激光雷达的距离，作为得分权重的参考
      std::vector<double> tmp_dis_set;

      for(int i = 0; i < segment_set.size(); ++i)
      {
        vector<Vector2<int>> grid_set;
        for(int j = 0; j < segment_set[i].size(); ++j)
        {
          const Vector2<double>& rPosition = segment_set[i][j].GetPosition();
          
          tmp_dis_set.push_back(rPosition.SquaredLength());

          // counterclockwise rotation and that rotation is about the origin (0, 0).
          Vector2<double> offset;
          offset.SetX(cosine * rPosition.GetX() - sine * rPosition.GetY());
          offset.SetY(-(sine * rPosition.GetX() + cosine * rPosition.GetY()));//转换到像素坐标系下，绕X轴旋转180度

          // have to compensate for the grid offset when getting the grid index
          Vector2<int> gridPoint = m_pGrid->WorldToGrid(offset + rGridOffset);

          grid_set.push_back(gridPoint);
        }

        std::sort(tmp_dis_set.begin(), tmp_dis_set.end());
        int dis_ind =  tmp_dis_set.size() / 2;
        double mid_dis = tmp_dis_set[dis_ind];

        dis_set.push_back(mid_dis);

        // double slope = cv_ws::cv_getLinePara(grid_set);
        // slope_set.push_back(slope);
        // std::cout << "slope " << i << ":" << slope << std::endl; 
      }

      //////////////////////////////////////////////////////
      // create lookup array for different angles
      double angle = 0.0;
      double startAngle = angleCenter - angleOffset;
      for (unsigned int angleIndex = 0; angleIndex < nAngles; angleIndex++)
      {
        angle = startAngle + angleIndex * angleResolution;
        ComputeOffsets(angleIndex, angle, localPoints_filter, pScan);
      }
      // assert(math::DoubleEqual(angle, angleCenter + angleOffset));
    }

  private:
    /**
     * Compute lookup value of points for given angle
     * @param angleIndex
     * @param angle
     * @param rLocalPoints
     */
    void ComputeOffsets(unsigned int angleIndex, double angle, const Pose2Vector& rLocalPoints, const LocalizedRangeScan& pScan)
    {
      m_ppLookupArray[angleIndex]->SetSize(static_cast<unsigned int>(rLocalPoints.size()));
      m_Angles.at(angleIndex) = angle;

      // set up point array by computing relative offsets to points readings
      // when rotated by given angle

      const Vector2<double>& rGridOffset = m_pGrid->GetCoordinateConverter()->GetOffset();

      double cosine = cos(angle);
      double sine = sin(angle);

      unsigned int readingIndex = 0;

      int* pAngleIndexPointer = m_ppLookupArray[angleIndex]->GetArrayPointer();

      const_forEach(Pose2Vector, &rLocalPoints)
      {
        const Vector2<double>& rPosition = iter->GetPosition();

        // counterclockwise rotation and that rotation is about the origin (0, 0).
        Vector2<double> offset;
        offset.SetX(cosine * rPosition.GetX() - sine * rPosition.GetY());
        offset.SetY(-(sine * rPosition.GetX() + cosine * rPosition.GetY()));//转换到像素坐标系下，绕X轴旋转180度

        // have to compensate for the grid offset when getting the grid index
        Vector2<int> gridPoint = m_pGrid->WorldToGrid(offset + rGridOffset);

        // use base GridIndex to ignore ROI
        int lookupIndex = m_pGrid->Grid<T>::GridIndex(gridPoint, false);

        pAngleIndexPointer[readingIndex] = lookupIndex;

        readingIndex++;
      }
      assert(readingIndex == rLocalPoints.size());
      std::cout << "readingIndex size = " << readingIndex << std::endl;
    }

    /**
     * Sets size of lookup table (resize if not big enough)
     * @param size
     */
    void SetSize(unsigned int size)
    {
      assert(size != 0);

      if (size > m_Capacity)
      {
        if (m_ppLookupArray != NULL)
        {
          DestroyArrays();
        }

        m_Capacity = size;
        m_ppLookupArray = new LookupArray*[m_Capacity];
        for (unsigned int i = 0; i < m_Capacity; i++)
        {
          m_ppLookupArray[i] = new LookupArray();
        }
      }

      m_Size = size;

      m_Angles.resize(size);
    }

    /**
     * Delete the arrays
     */
    void DestroyArrays()
    {
      for (unsigned int i = 0; i < m_Capacity; i++)
      {
        delete m_ppLookupArray[i];
      }

      delete[] m_ppLookupArray;
      m_ppLookupArray = NULL;
    }

  private:
    Grid<T>* m_pGrid;

    unsigned int m_Capacity;
    unsigned int m_Size;

    LookupArray **m_ppLookupArray;

    // for sanity check
    std::vector<double> m_Angles;

    std::vector<unsigned int> seg_counter;

    std::vector<double> slope_set;
    std::vector<double> dis_set;

    float m_dis_threshold;
    int m_batch_size;
  };  // class GridIndexLookup

  class CorrelationGrid : public Grid<uint8_t>
  {
  public:
    /**
     * Destructor
     */
    virtual ~CorrelationGrid()
    {
       if(m_pKernel)
         delete [] m_pKernel;
      //  std::cout<<"Deconstruction CorrelationGrid"<<std::endl;
    }

  public:
    /**
     * Create a correlation grid of given size and parameters
     * @param width
     * @param height
     * @param resolution
     * @param smearDeviation
     * @return correlation grid
     */
    static CorrelationGrid* CreateGrid(int width,
                                       int height,
                                       double resolution,
                                       double smearDeviation)
    {
      assert(resolution != 0.0);

      // +1 in case of roundoff
      // unsigned int borderSize = GetHalfKernelSize(smearDeviation, resolution) + 1;

      unsigned int borderSize = 0;

      CorrelationGrid* pGrid = new CorrelationGrid(width, height, borderSize, resolution, smearDeviation);

      return pGrid;
    }

    /**
     * Gets the index into the data pointer of the given grid coordinate
     * @param rGrid
     * @param boundaryCheck
     * @return grid index
     */
    virtual int GridIndex(const Vector2<int>& rGrid, bool boundaryCheck = false) const
    {
      int x = rGrid.GetX();
      int y = rGrid.GetY();

      return Grid<uint8_t>::GridIndex(Vector2<int>(x, y), boundaryCheck);
    }

    /**
     * Get the Region Of Interest (ROI)
     * @return region of interest
     */
    inline const Rectangle2<int>& GetROI() const
    {
      return m_Roi;
    }

    /**
     * Sets the Region Of Interest (ROI)
     * @param roi
     */
    inline void SetROI(const Rectangle2<int>& roi)
    {
      m_Roi = roi;
    }

    /**
     * Smear cell if the cell at the given point is marked as "occupied"
     * @param rGridPoint
     */
    inline void SmearPoint(const int& index)
    {
      assert(m_pKernel != NULL);

      if (GetDataPointer()[index] != GridStates_Occupied)
      {
        return;
      }

      int halfKernel = m_KernelSize / 2;

      // apply kernel
      for (int j = -halfKernel; j <= halfKernel; j++)
      {
        int offset = index + j * GetWidth();

        if(offset < halfKernel || offset > GetWidth()*GetHeight()-halfKernel)
        {
          continue;
        }

        uint8_t* pGridAdr = GetDataPointer() + offset;

        int kernelConstant = (halfKernel) + m_KernelSize * (j + halfKernel);

        // if a point is on the edge of the grid, there is no problem
        // with running over the edge of allowable memory, because
        // the grid has margins to compensate for the kernel size
        for (int i = -halfKernel; i <= halfKernel; i++)
        {
          int kernelArrayIndex = i + kernelConstant;

          uint8_t kernelValue = m_pKernel[kernelArrayIndex];
          if (kernelValue > pGridAdr[i])
          {
            // kernel value is greater, so set it to kernel value
            pGridAdr[i] = kernelValue;
          }
        }
      }
    }

  protected:
    /**
     * Constructs a correlation grid of given size and parameters
     * @param width
     * @param height
     * @param borderSize
     * @param resolution
     * @param smearDeviation
     */
    /*               Grid构造示意图
    ---------------------------------------------
    *                                            *
    *                                            *
    *                                            *
    *                                            *
    *                      width                 *
    *               -----------------            *
    *               *               *            *
    *               *               * height     *
    *               *               *            *
    *  boardersize  *               *            *
    * <------------>-----------------            *
    *                                            *
    *                                            *
    *                                            *
    *                                            *
    *--------------------------------------------*
     */
    CorrelationGrid(unsigned int width, unsigned int height, unsigned int borderSize,
                    double resolution, double smearDeviation)
      : Grid<uint8_t>(width + borderSize * 2, height + borderSize * 2)
      , m_SmearDeviation(smearDeviation)
    {
      GetCoordinateConverter()->SetScale(1.0 / resolution);

      // setup region of interest
       m_Roi = Rectangle2<int>(borderSize, borderSize, width, height);

      // calculate kernel
       CalculateKernel();
    }

    /**
     * Sets up the kernel for grid smearing.
     */
    virtual void CalculateKernel()
    {
      double resolution = GetResolution();

      assert(resolution != 0.0);
      assert(m_SmearDeviation != 0.0);

      // min and max distance deviation for smearing;
      // will smear for two standard deviations, so deviation must be at least 1/2 of the resolution
      const double MIN_SMEAR_DISTANCE_DEVIATION = 0.5 * resolution;
      const double MAX_SMEAR_DISTANCE_DEVIATION = 10 * resolution;

      // check if given value too small or too big
      if (!math::InRange(m_SmearDeviation, MIN_SMEAR_DISTANCE_DEVIATION, MAX_SMEAR_DISTANCE_DEVIATION))
      {
        std::stringstream error;
        error << "Mapper Error:  Smear deviation too small:  Must be between "
              << MIN_SMEAR_DISTANCE_DEVIATION
              << " and "
              << MAX_SMEAR_DISTANCE_DEVIATION;
        throw std::runtime_error(error.str());
      }

      // NOTE:  Currently assumes a two-dimensional kernel

      // +1 for center
      m_KernelSize = 2 * GetHalfKernelSize(m_SmearDeviation, resolution) + 1;

      // allocate kernel
      m_pKernel = new uint8_t[m_KernelSize * m_KernelSize];
      if (m_pKernel == NULL)
      {
        throw std::runtime_error("Unable to allocate memory for kernel!");
      }

      // calculate kernel
      int halfKernel = m_KernelSize / 2;
      for (int i = -halfKernel; i <= halfKernel; i++)
      {
        for (int j = -halfKernel; j <= halfKernel; j++)
        {
#ifdef WIN32
          double distanceFromMean = _hypot(i * resolution, j * resolution);
#else
          double distanceFromMean = hypot(i * resolution, j * resolution);
#endif
          double z = exp(-0.5 * pow(distanceFromMean / m_SmearDeviation, 2));

          unsigned int kernelValue = static_cast<unsigned int>(math::Round(z * GridStates_Occupied));
          assert(math::IsUpTo(kernelValue, static_cast<unsigned int>(255)));

          int kernelArrayIndex = (i + halfKernel) + m_KernelSize * (j + halfKernel);
          m_pKernel[kernelArrayIndex] = static_cast<uint8_t>(kernelValue);
        }
      }
    }

    /**
     * Computes the kernel half-size based on the smear distance and the grid resolution.
     * Computes to two standard deviations to get 95% region and to reduce aliasing.
     * @param smearDeviation
     * @param resolution
     * @return kernel half-size based on the parameters
     */
    static int GetHalfKernelSize(double smearDeviation, double resolution)
    {
      assert(resolution != 0.0);

      return static_cast<int>(math::Round(2.0 * smearDeviation / resolution));
    }

  private:
    /**
     * The point readings are smeared by this value in X and Y to create a smoother response.
     * Default value is 0.03 meters.
     */
    double m_SmearDeviation;

    // Size of one side of the kernel
     int m_KernelSize;

    // Cached kernel for smearing
     uint8_t* m_pKernel;

    // region of interest
     Rectangle2<int> m_Roi;
  };  // CorrelationGrid

}//namespace karto

#endif