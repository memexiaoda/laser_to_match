#ifndef PARAMETER_H
#define PARAMETER_H

#include "Types.h"

namespace karto
{

template <class T>
class Parameter
{
public:
  Parameter(){}

  Parameter(T val)
  {
  	SetValue(val);
  }

  void SetValue(T val)
  {
	m_data = val;
  }
	
  T GetValue()
  {
	return m_data;
  }
	
private:
  T m_data;

};

class Mapper
{
	public:
    Mapper():m_pCorrelationSearchSpaceDimension(NULL),
             m_pCorrelationSearchSpaceResolution(NULL),
             m_pCorrelationSearchSpaceSmearDeviation(NULL),
             m_pDistanceVariancePenalty(NULL),
             m_pAngleVariancePenalty(NULL),
             m_pFineSearchAngleOffset(NULL),
             m_pCoarseSearchAngleOffset(NULL),
             m_pMinimumAnglePenalty(NULL),
             m_pMinimumDistancePenalty(NULL),
             m_pUseResponseExpansion(NULL){}




	Parameter<double>* m_pCorrelationSearchSpaceDimension;
	Parameter<double>* m_pCorrelationSearchSpaceResolution;
    Parameter<double>* m_pCorrelationSearchSpaceSmearDeviation;

	Parameter<double>* m_pDistanceVariancePenalty;
	Parameter<double>* m_pAngleVariancePenalty;
	Parameter<double>* m_pFineSearchAngleOffset;
	Parameter<double>* m_pCoarseSearchAngleOffset;
	Parameter<double>* m_pCoarseAngleResolution;
	Parameter<double>* m_pMinimumAnglePenalty;
	Parameter<double>* m_pMinimumDistancePenalty;
	Parameter<bool>* m_pUseResponseExpansion;
	
	public:
	~Mapper()
	{
		if(m_pCorrelationSearchSpaceDimension)
			delete m_pCorrelationSearchSpaceDimension;
		if(m_pCorrelationSearchSpaceResolution)
			delete m_pCorrelationSearchSpaceResolution;
		if(m_pCorrelationSearchSpaceSmearDeviation)
			delete m_pCorrelationSearchSpaceSmearDeviation;

		if(m_pCoarseAngleResolution)
			delete m_pCoarseAngleResolution;
		if(m_pDistanceVariancePenalty)
			delete m_pDistanceVariancePenalty;
		if(m_pAngleVariancePenalty)
			delete m_pAngleVariancePenalty;
		if(m_pFineSearchAngleOffset)
			delete m_pFineSearchAngleOffset;
		if(m_pCoarseSearchAngleOffset)
			delete m_pCoarseSearchAngleOffset;
		if(m_pCoarseAngleResolution)
			delete m_pCoarseAngleResolution;
		if(m_pMinimumAnglePenalty)
			delete m_pMinimumAnglePenalty;
		if(m_pMinimumDistancePenalty)
			delete m_pMinimumDistancePenalty;
		if(m_pUseResponseExpansion)
			delete m_pUseResponseExpansion;
	}
	
	
	public:
	
	// Setting Scan Matcher Parameters from the Parameter Server

	void setParamCorrelationSearchSpaceDimension(double val)
	{
		m_pCorrelationSearchSpaceDimension = new Parameter<double>();
		m_pCorrelationSearchSpaceDimension->SetValue(val);
	}

	void setParamCorrelationSearchSpaceResolution(double val)
	{
		m_pCorrelationSearchSpaceResolution = new Parameter<double>();
		m_pCorrelationSearchSpaceResolution->SetValue(val);
	}

	void setParamCorrelationSearchSpaceSmearDeviation(double val)
	{
		m_pCorrelationSearchSpaceSmearDeviation = new Parameter<double>();
		m_pCorrelationSearchSpaceSmearDeviation->SetValue(val);
	}


    void setParamDistanceVariancePenalty(double val)
	{
		m_pDistanceVariancePenalty = new Parameter<double>();
		m_pDistanceVariancePenalty->SetValue(val);
	}


    void setParamAngleVariancePenalty(double val)
	{
		m_pAngleVariancePenalty = new Parameter<double>();
		m_pAngleVariancePenalty->SetValue(val);
	}


    void setParamFineSearchAngleOffset(double val)
	{
		m_pFineSearchAngleOffset = new Parameter<double>();
		m_pFineSearchAngleOffset->SetValue(val);
	}


    void setParamCoarseSearchAngleOffset(double val)
	{
		m_pCoarseSearchAngleOffset = new Parameter<double>();
		m_pCoarseSearchAngleOffset->SetValue(val);
	}


    void setParamCoarseAngleResolution(double val)
	{
		m_pCoarseAngleResolution = new Parameter<double>();
		m_pCoarseAngleResolution->SetValue(val);
	}


    void setParamMinimumAnglePenalty(double val)
	{
		m_pMinimumAnglePenalty = new Parameter<double>();
		m_pMinimumAnglePenalty->SetValue(val);
	}


    void setParamMinimumDistancePenalty(double val)
	{
		m_pMinimumDistancePenalty = new Parameter<double>();
		m_pMinimumDistancePenalty->SetValue(val);
	}


    void setParamUseResponseExpansion(bool val)
	{
		m_pUseResponseExpansion = new Parameter<bool>();
		m_pUseResponseExpansion->SetValue(val);
	}
};
}//namespace karto

#endif
