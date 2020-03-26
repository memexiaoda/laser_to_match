#include "ros/ros.h"
#include "ros/console.h"

#include <signal.h>

#include <boost/circular_buffer.hpp>
#include "tf/transform_broadcaster.h"
#include "tf/transform_listener.h"
#include "tf/message_filter.h"
#include "tf/tf.h"
#include "message_filters/subscriber.h"

#include "nav_msgs/MapMetaData.h"
#include "std_msgs/Float32.h"
#include "sensor_msgs/LaserScan.h"
#include "nav_msgs/GetMap.h"
#include "nav_msgs/SetMap.h"
#include "std_srvs/Empty.h"

#include "tf2/LinearMath/Transform.h"
#include "tf2/convert.h"
// #include "tf2/utils.h"
//#include "tf2_geometry_msgs/tf2_geometry_msgs.h"
#include "tf2_ros/buffer.h"
#include "tf2_ros/message_filter.h"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_ros/transform_listener.h"

#include "location_match_score/scan_match.h"
#include "location_match_score/Pose.h"
#include "location_match_score/Grid.h"

// #include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>

#include <string>
#include <map>
#include <vector>
#include <memory>

#include <eigen3/Eigen/Dense>
#include "location_match_score/segment.h"

// compute linear index for given map coords
// #define MAP_IDX(sx, i, j) ((sx) * (j) + (i))

typedef enum
{
  GridStates_Unknown = 0,
  GridStates_Occupied = 100,
  GridStates_Free = 0
} GridStates;

void sigintHandler(int sig)
{
  // Save latest pose as we're shutting down.
  // amcl_node_ptr->savePoseToServer();
  ros::shutdown();
}

class LocationMatchScore
{
  public:
    LocationMatchScore();
    ~LocationMatchScore();

    void laserCallback(const sensor_msgs::LaserScan::ConstPtr& scan);
    void requestMap();

  private:
    bool getLaserPose(karto::Pose2& karto_pose, const ros::Time& t, std::string frame_id);
    std::shared_ptr<karto::LaserRangeFinder> getLaser(const sensor_msgs::LaserScan::ConstPtr& scan);
    bool addScan(const sensor_msgs::LaserScan::ConstPtr& scan, karto::Pose2& karto_pose);
    bool updateMap();

    void mapReceived(const nav_msgs::OccupancyGridConstPtr& msg);
    void handleMapMessage(const nav_msgs::OccupancyGrid& msg);
    void convertMap( const nav_msgs::OccupancyGrid& map_msg );

    // ROS handles
    ros::NodeHandle nh_, private_nh_;
    ros::Publisher response_pub_;
    struct TransformListenerWrapper : public tf::TransformListener
    {
      // inline tf2_ros::Buffer &getBuffer() {return tf2_buffer_;}
    };

    TransformListenerWrapper* tf_;

    message_filters::Subscriber<sensor_msgs::LaserScan>* laser_scan_sub_;
    tf::MessageFilter<sensor_msgs::LaserScan>* laser_scan_filter_;

    ros::Subscriber map_sub_;

    // The map that will be published / send to service callers
    nav_msgs::GetMap::Response map_;
    bool first_map_only_;

    // Storage for ROS parameters
    std::string odom_frame_;
    std::string map_frame_;
    std::string base_frame_;
    std::string laser_scan_frame_id;
    int throttle_scans_;
    ros::Duration map_update_interval_;
    double resolution_;
    boost::mutex map_mutex_;
    boost::mutex map_to_odom_mutex_;

    // Karto bookkeeping
    karto::Mapper* mapper_;
    karto::ScanMatcher* scanmatcher_;
    karto::CorrelationGrid* m_pCorrelationGrid_;

    std::map<std::string, std::shared_ptr<karto::LaserRangeFinder> > lasers_;

    // Internal state
    bool got_map_;
    int laser_count_;

    double min_threshold_laser_;
    double max_threshold_laser_;

    double rangeThreshold_;

    bool first_map_received_;
    bool use_map_topic_;
    tf::Transform map_to_odom_;

    boost::recursive_mutex configuration_mutex_;

    double response;
    boost::circular_buffer<double> cir_buf;
};

boost::scoped_ptr<LocationMatchScore> LocationMatchScore_node_ptr;

LocationMatchScore::LocationMatchScore() :
        got_map_(false),
        laser_count_(0),
        private_nh_("~"),
        mapper_(NULL),
        scanmatcher_(NULL),
        m_pCorrelationGrid_(NULL)
{
  boost::recursive_mutex::scoped_lock l(configuration_mutex_);
  map_to_odom_.setIdentity();
  // Retrieve parameters
  if(!private_nh_.getParam("odom_frame", odom_frame_))
    odom_frame_ = "odom";
  if(!private_nh_.getParam("map_frame", map_frame_))
    map_frame_ = "map";
  if(!private_nh_.getParam("base_frame", base_frame_))
    base_frame_ = "base_link";
  if(!private_nh_.getParam("throttle_scans", throttle_scans_))
    throttle_scans_ = 5;

  if(!private_nh_.getParam("min_threshold_laser", min_threshold_laser_))
    min_threshold_laser_ = 0.2;

  if(!private_nh_.getParam("max_threshold_laser", max_threshold_laser_))
    max_threshold_laser_ = 20.0;

  private_nh_.param("first_map_only", first_map_only_, true);
  private_nh_.param("use_map_topic", use_map_topic_, false);
  private_nh_.param("first_map_received", first_map_received_, false);

  tf_ = new TransformListenerWrapper();

  // Initialize Karto structures
  mapper_ = new karto::Mapper();

  double correlation_search_space_dimension_;
  private_nh_.param("correlation_search_space_dimension", correlation_search_space_dimension_, 0.04);
  mapper_->setParamCorrelationSearchSpaceDimension(correlation_search_space_dimension_);

  double correlation_search_space_resolution_;
  private_nh_.param("correlation_search_space_resolution", correlation_search_space_resolution_, 0.01);
  mapper_->setParamCorrelationSearchSpaceResolution(correlation_search_space_resolution_);

  double correlation_search_space_smear_deviation_;
  private_nh_.param("correlation_search_space_smear_deviation", correlation_search_space_smear_deviation_, 0.04);
  mapper_->setParamCorrelationSearchSpaceSmearDeviation(correlation_search_space_smear_deviation_);

  // Setting Correlation Parameters, Loop Closure Parameters from the Parameter Server

  // Setting Scan Matcher Parameters from the Parameter Server

  double distance_variance_penalty_;
  private_nh_.param("distance_variance_penalty", distance_variance_penalty_, 0.3);
  mapper_->setParamDistanceVariancePenalty(distance_variance_penalty_);

  double angle_variance_penalty_;
  private_nh_.param("angle_variance_penalty", angle_variance_penalty_, 0.349);
  mapper_->setParamAngleVariancePenalty(angle_variance_penalty_);

  double fine_search_angle_offset_;
  private_nh_.param("fine_search_angle_offset", fine_search_angle_offset_, 0.00349);
  mapper_->setParamFineSearchAngleOffset(fine_search_angle_offset_);

  double coarse_search_angle_offset_;
  private_nh_.param("coarse_search_angle_offset", coarse_search_angle_offset_, 0.08725);//0.1745
  mapper_->setParamCoarseSearchAngleOffset(coarse_search_angle_offset_);

  double coarse_angle_resolution_;
  private_nh_.param("coarse_angle_resolution", coarse_angle_resolution_, 0.0349);
  mapper_->setParamCoarseAngleResolution(coarse_angle_resolution_);

  double minimum_angle_penalty_;
  private_nh_.param("minimum_angle_penalty", minimum_angle_penalty_, 0.9);
  mapper_->setParamMinimumAnglePenalty(minimum_angle_penalty_);

  double minimum_distance_penalty_;
  private_nh_.param("minimum_distance_penalty", minimum_distance_penalty_, 0.5);
  mapper_->setParamMinimumDistancePenalty(minimum_distance_penalty_);

  bool use_response_expansion;
  private_nh_.param("use_response_expansion", use_response_expansion, true);
  mapper_->setParamUseResponseExpansion(use_response_expansion);

  private_nh_.param("rangeThreshold", rangeThreshold_, 15.0);

  private_nh_.param("use_map_topic", use_map_topic_, true);

  response_pub_ = nh_.advertise<std_msgs::Float32>("match_response",1,true);


  laser_scan_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(nh_, "scan", 5);

  laser_scan_filter_ = new tf::MessageFilter<sensor_msgs::LaserScan>(*laser_scan_sub_, 
                                                                      *tf_, 
                                                                      odom_frame_, 
                                                                      5);

  laser_scan_filter_->registerCallback(boost::bind(&LocationMatchScore::laserCallback, this, _1));

  if(use_map_topic_)
  {
    map_sub_ = nh_.subscribe("map", 1, &LocationMatchScore::mapReceived, this);
    ROS_INFO("Subscribed to map topic.");
  }else{
    requestMap();
  }
  cir_buf.resize(5);

}

LocationMatchScore::~LocationMatchScore()
{

  if (laser_scan_sub_)
    delete laser_scan_sub_;

  if (laser_scan_filter_)
    delete laser_scan_filter_;

  if(scanmatcher_)
    delete scanmatcher_;

  if (mapper_)
    delete mapper_;
  if(m_pCorrelationGrid_)
    delete m_pCorrelationGrid_;

  // if(lasers_[laser_scan_frame_id])
  //   delete lasers_[laser_scan_frame_id];

  delete tf_;
  // TODO: delete the pointers in the lasers_ map; not sure whether or not
  // I'm supposed to do that.
}


std::shared_ptr<karto::LaserRangeFinder>
LocationMatchScore::getLaser(const sensor_msgs::LaserScan::ConstPtr& scan)
{

  if(lasers_.find(scan->header.frame_id) == lasers_.end())
  {    
    std::string name = scan->header.frame_id;
    std::shared_ptr<karto::LaserRangeFinder> laser = 
      karto::LaserRangeFinder::CreateLaserRangeFinder();

    laser->SetMinimumRange(scan->range_min);
    laser->SetMaximumRange(scan->range_max);
    laser->SetMinimumAngle(scan->angle_min);
    laser->SetMaximumAngle(scan->angle_max);
    laser->SetAngularResolution(scan->angle_increment);

    lasers_[scan->header.frame_id] = laser;
  }

  return lasers_[scan->header.frame_id];
  
}


bool
LocationMatchScore::getLaserPose(karto::Pose2& karto_pose, const ros::Time& t, std::string frame_id)
{
  // Get the laser's pose

  tf::Stamped<tf::Pose> ident (tf::Transform(tf::createQuaternionFromRPY(0,0,0),
                                           tf::Vector3(0,0,0)), t, frame_id);

  tf::Stamped<tf::Transform> laser_pose;
  try
  {
     tf_->transformPose(map_frame_, ident, laser_pose);
  }
  catch(tf2::TransformException e)
  {
    ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
    return false;
  }
  double yaw = tf::getYaw(laser_pose.getRotation());

  karto_pose = karto::Pose2(laser_pose.getOrigin().x(),
                           laser_pose.getOrigin().y(),
                           yaw);
  // ROS_INFO("laser pose: x = %f, y = %f, yaw = %f ", 
  //          laser_pose.getOrigin().x(),
  //          laser_pose.getOrigin().y(),
  //          yaw);

  return true;
}


void
LocationMatchScore::laserCallback(const sensor_msgs::LaserScan::ConstPtr& scan)
{
  if(!first_map_received_)
    return;

  boost::recursive_mutex::scoped_lock lc(configuration_mutex_);

  laser_scan_frame_id = scan->header.frame_id;
  //ROS_INFO("laser_scan_frame_id, %s", laser_scan_frame_id.c_str());
  // Check whether we know about this laser yet
  // std::shared_ptr<karto::LaserRangeFinder> laser = getLaser(scan);

  // if(!laser)
  // {
  //   ROS_WARN("Failed to create laser device for %s; discarding scan",
  //      scan->header.frame_id.c_str());
  //   return;
  // }

  if(scanmatcher_)
  {
    laser_count_++;
    if ((laser_count_ % throttle_scans_) != 0)
    return;

    // static ros::Time last_map_update(0,0);

    karto::Pose2 laser_pose;
    if(addScan(scan, laser_pose))
    {
      cir_buf.push_back(response);
      double tmp_sum = 0;
      for(int i = 0; i < cir_buf.size(); ++i)
      {
         tmp_sum += cir_buf[i];
      }
      double avr_res = tmp_sum / cir_buf.size();
      std_msgs::Float32 res;
      res.data = avr_res;
      response_pub_.publish(res);
      ROS_DEBUG("added scan at pose: %.3f %.3f %.3f", laser_pose.GetX(),laser_pose.GetY(),laser_pose.GetHeading());
    }
  }
  
}


bool
LocationMatchScore::addScan(const sensor_msgs::LaserScan::ConstPtr& scan, 
                   karto::Pose2& karto_pose)
{
  if(!getLaserPose(karto_pose, scan->header.stamp, scan->header.frame_id))
     return false;
  
  // Create a vector of doubles for karto
  std::vector<Eigen::Vector2f> reading_set;

  int index = 0;

  for(std::vector<float>::const_iterator it = scan->ranges.begin();
    it != scan->ranges.end();
      ++it, ++index)
  {
    if(*it > min_threshold_laser_ && *it < max_threshold_laser_)
    {
      float angle = scan->angle_min + index * scan->angle_increment;

      float x = *it * cos(angle);
      float y = *it * sin(angle);

      Eigen::Vector2f val(x, y);
      reading_set.push_back(val);  
    }else{
      continue;
    }
    
  }

  std::cout << "reading_set size = " << reading_set.size() << std::endl;
  // ROS_INFO("number of laser points = %d", count);
  // ROS_INFO("laser readings copied is done!");

  // std::vector<std::vector<Eigen::Vector2f> > seg_set;
  // seg_set = segment(reading_set);
   ros::Time begin = ros::Time::now();

  karto::LocalizedRangeScan range_scan;
  range_scan.SetRangeReadings(reading_set);
  range_scan.SetLaserRangeFinder(lasers_[scan->header.frame_id]);
  range_scan.SetOdometricPose(karto_pose);
  range_scan.SetCorrectedPose(karto_pose);
  karto::Pose2 rMean;
  karto::Matrix3 rCovariance;
  response = scanmatcher_->MatchScan(range_scan, rMean, rCovariance, false, false);

  ROS_INFO("score %f", response);

  ros::Time end = ros::Time::now();
  ros::Duration calcu_time = end - begin;

  ROS_INFO("Calculate time = %f", calcu_time.toSec());

  return true;
}

void
LocationMatchScore::mapReceived(const nav_msgs::OccupancyGridConstPtr& msg)
{
  if( first_map_only_ && first_map_received_ ) {
    return;
  }
  ROS_INFO("map received");
  handleMapMessage( *msg );

  ROS_INFO("map has been handled!");


  first_map_received_ = true;
}

void
LocationMatchScore::requestMap()
{
  boost::recursive_mutex::scoped_lock ml(configuration_mutex_);

  // get map via RPC
  nav_msgs::GetMap::Request  req;
  nav_msgs::GetMap::Response resp;
  ROS_INFO("Requesting the map...");
  while(!ros::service::call("static_map", req, resp))
  {
    ROS_WARN("Request for map failed; trying again...");
    ros::Duration d(0.5);
    d.sleep();
  }
  handleMapMessage( resp.map );

  first_map_received_ = true;
}

void
LocationMatchScore::handleMapMessage(const nav_msgs::OccupancyGrid& msg)
{
  boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);

  convertMap(msg);
}

void
LocationMatchScore::convertMap(const nav_msgs::OccupancyGrid& map_msg)
{
  ROS_INFO("Mapper SmearDeviation: %f", mapper_->m_pCorrelationSearchSpaceSmearDeviation->GetValue());
  m_pCorrelationGrid_ = karto::CorrelationGrid::CreateGrid(map_msg.info.width, map_msg.info.height, map_msg.info.resolution, mapper_->m_pCorrelationSearchSpaceSmearDeviation->GetValue());

  ROS_INFO("start parse map message!");

  ROS_INFO("map width = %d", map_msg.info.width);
  ROS_INFO("map height = %d", map_msg.info.height);

  // Occupancy state (-1 = free, 0 = unknown, +1 = occ) from amcl
  int count_of_occupied = 0;
  int _index = 0;
  for(int j = map_msg.info.height - 1; j >= 0; j--)
  {
    for(int i = 0;i< map_msg.info.width;i++)
    {
      int tt = j * map_msg.info.width + i;
      if(map_msg.data[tt] == 0)
      {
        m_pCorrelationGrid_->GetDataPointer()[_index]  = GridStates_Free;//free
      }
      else if(map_msg.data[tt] == 100)
      {
        m_pCorrelationGrid_->GetDataPointer()[_index] = GridStates_Occupied;//occ
        count_of_occupied++;
      }
      else{
        m_pCorrelationGrid_->GetDataPointer()[_index]= GridStates_Unknown;//unknown
      }
      _index++; 
  }
  }
  std::cout << "count_of_occupied = " << count_of_occupied << std::endl;

  karto::Vector2<double> offset;
  offset.SetX(map_msg.info.origin.position.x);
  offset.SetY(map_msg.info.origin.position.y);

  std::cout << "offset.SetX " <<  offset.GetX();
  std::cout << "offset.SetY " << offset.GetY();

  m_pCorrelationGrid_->GetCoordinateConverter()->SetOffset(offset);

  std::cout << "m_pCorrelationGrid_ "<< m_pCorrelationGrid_<<std::endl;

  if(m_pCorrelationGrid_)
  {
    // std::cout <<"create scanmatcher!!!"<<std::endl;
    scanmatcher_ = karto::ScanMatcher::Create(mapper_,
                                             mapper_->m_pCorrelationSearchSpaceDimension->GetValue(),
                                             mapper_->m_pCorrelationSearchSpaceResolution->GetValue(),
                                             mapper_->m_pCorrelationSearchSpaceSmearDeviation->GetValue(),
                                             rangeThreshold_,
                                             m_pCorrelationGrid_);
  }

  ROS_INFO("Grid is ok to use!!!"); 
}


int
main(int argc, char** argv)
{
  ros::init(argc, argv, "location_match_score", ros::init_options::NoSigintHandler);
  // ros::NodeHandle nh;

  // signal(SIGINT, sigintHandler);

  // LocationMatchScore_node_ptr.reset(new LocationMatchScore());

  LocationMatchScore lm;

  ros::spin();

  // LocationMatchScore_node_ptr.reset();

  return 0;
}
