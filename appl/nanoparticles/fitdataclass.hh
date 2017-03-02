// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARAMETERESTIMATION_FITDATA_HH
#define DUNE_DYCAP_PARAMETERESTIMATION_FITDATA_HH


#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <map>

#include <sys/types.h>
#include <sys/stat.h>

#include <src/estimation/utilities.hh>
#include <src/estimation/fitclass.hh>

/** \brief FitDataClass for parameter estimation
 */
class DycapFitDataClass : public FitDataClass
{


public:
  DycapFitDataClass(string infilename, int verbosity_ = 0) :
    modelFileName_("model.dat"), outputCounter_(0), verbosity(verbosity_)
  {
    ReadData(infilename);
  }

  DycapFitDataClass(int verbosity_ = 0) :
    modelFileName_("model.dat"), outputCounter_(0), verbosity(verbosity_)
  {
  }


  //! generates solution in each iteration step
  void OutputSolution(const vector<FLOAT> yModel)
  {
    outputCounter_++;

    std::string solution_path = "model_output";
    mkdir(solution_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::stringstream fileName;

    std::ostringstream s;
    s << "./" + solution_path + "/" << outputCounter_  << "_" << modelFileName_;
    srand((unsigned)time(NULL));
    ofstream outfile;
    outfile.open(s.str().c_str());
    if (!outfile.is_open())
      {
        cerr << "Not able to open file " << fileName.str().c_str() << "!" << endl;
        throw error_abort();
      }

    outfile.precision(14);
    for(DataMapType::iterator element=measurements_.begin();element!=measurements_.end();++element)
      {
        outfile << element->first;
        for(size_t j=0;j<element->second.size();++j)
          {
            if (element->second[j].active)
#ifdef GENERATE_VALUES
              if (element->second[j].stdDev<0)
                outfile << "\t" << yModel[element->second[j].index] * (1.-element->second[j].stdDev*((double(rand())/RAND_MAX)-0.5));
              else
                outfile << "\t" << yModel[element->second[j].index] + element->second[j].stdDev*((double(rand())/RAND_MAX)-0.5);
#else
            outfile << "\t" << yModel[element->second[j].index];
#endif
            else
              outfile << "\t" << -1e+300;
            outfile << "\t" << element->second[j].value;
          }
        outfile << endl;
      }

    outfile.close();
  }



  INT GetNumValues()
  {
    return(data_.size());
  }

  //! get vector with all timesteps
  vector<FLOAT> & GetTimes()
  {
    return(time_);
  }


  FLOAT GetValue(INT i)
  {
    return(data_[i].value);
  }

  //! get measurements with concrete time
  const vector<measurementType> GetMeasurements(FLOAT time)
  {
    DataMapType::iterator measurement=measurements_.find(time);
    if (measurement!=measurements_.end())
      return(measurement->second);
    else
      {
        vector<measurementType> fm = (measurements_.begin())->second;
        for(size_t i=0;i<fm.size();++i)
          fm[i].value = -5;
        return fm;
      }
  }

  //! number of all experiments
  int getNumberOfExperiments()
  {
    return numberofexperiments;
  }

  //! get measurements between two times
  const DataMapType GetMeasurementArray(FLOAT oldTime,FLOAT newTime)
  {
    DataMapType::iterator lower=measurements_.upper_bound(oldTime);
    DataMapType::iterator upper=measurements_.lower_bound(newTime);
    if ((lower==measurements_.end()) || (lower->first > newTime))
      return(DataMapType());
    if ((upper!=measurements_.end()) && ((upper==lower) || (fabs(upper->first-newTime)<1e-30)))
      ++upper;

    return(DataMapType(lower,upper));
  }

  FLOAT GetStdDev(INT i)
  {
    return(data_[i].stdDev);
  }

  string GetType(INT i)
  {
    return(data_[i].type);
  }

  void ReadData(string infilename)
  {
    // read measured data field
    ifstream infile(infilename.c_str());
    if (!infile)
      {
        cerr << "DycapFitDataClass::ReadData Error: Could not open data file " << infilename << "!" << endl;
        throw error_abort();
      }
    else
      {
        if (verbosity)
          std::cout << "File " << infilename.c_str() << " was open" << std::endl;
      }

    SkipComment(infile);
    // first line is number of measurements (timesteps)!
    FLOAT numData;
    infile >> numData;
    if (verbosity)
      std::cout << "Number of measurements (timesteps)" << numData << std::endl;

    SkipComment(infile);
    // second line is number of experiments!
    infile >> numberofexperiments;

    if (verbosity)
      std::cout << "Number of experiments is " << numberofexperiments << std::endl;

    char buffer[256];
    SkipComment(infile);
    infile.getline(buffer,256);
    istringstream instream(buffer);
    // type of variables
    string type;
    vector<string> types_;
    while (instream>>type)
      {
        for (size_t i=0;i<type.size();++i)
          type[i]=tolower(type[i]);
        types_.push_back(type);
      }

    if (verbosity){
      std::cout << "variable types ";
      for (auto it = types_.begin(); it != types_.end(); ++it )
        std::cout << *it << " ";
      std::cout << "\n";
    }

    SkipComment(infile);
    infile.getline(buffer,256);
    instream.clear();
    instream.str(buffer);
    FLOAT value;
    vector<FLOAT> stdDev;
    while (instream>>value)
      {
        char relative;
        if (instream.get(relative))
          {
            if (tolower(relative)=='r')
              value=-value;
            else
              instream.unget();
          }
        stdDev.push_back(value);
      }
    if (verbosity){
      std::cout << "standard deviation ";
      for (auto it = stdDev.begin(); it != stdDev.end(); ++it )
        std::cout << *it << " ";
      std::cout << "\n\n";
    }


    std::string fname = infilename;
    std::ofstream soutput;
    fname.resize(fname.length()-4);
    fname.append("_data.dat");
    soutput.open(fname.c_str());
    soutput << std::setprecision(14) << std::scientific;

    measurements_.clear();
    FLOAT time=0.;
    for(INT line=0;line<numData;++line)
      {

        SkipComment(infile);
        infile.getline(buffer,256);
        instream.clear();

        instream.str(buffer);
        FLOAT lastTime=time;

        instream >> time;
        soutput << time << " ";

        vector<measurementType> measurements;

        int lindex(0);
        while (instream>>value)
          {
            measurementType measurement;

            // be carefull, negative values means inactive!!
            if ( (fabs(stdDev[measurements.size()])>1e-50) && value>=0)
              measurement.active=true;
            else
              measurement.active=false;

            measurement.type=types_[measurements.size()];
            measurement.localindex=lindex%numberofexperiments;

            measurement.value = value;

            if (value>=0.)
              soutput << measurement.value << " ";
            else
              soutput << 0.0 << " ";

            if (stdDev[measurements.size()]>0.)
              measurement.stdDev=stdDev[measurements.size()];
            else
              measurement.stdDev=measurement.value*-stdDev[measurements.size()];
            measurements.push_back(measurement);
            lindex++;
          }


        measurements_.insert(make_pair(time,measurements));

        if (verbosity)
          {
            std::cout <<  right << fixed << setw(8) << setprecision(0) << time << "    ";
            auto mv = GetMeasurements(time);
            for (auto it = mv.begin(); it != mv.end(); ++it )
              std::cout <<  left << fixed << setw(12) << setprecision(3) << (*it).value;
            std::cout << "\n";
          }

        soutput << "\n";
      }
    std::cout << "\n";

    soutput.close();
    infile.close();
    time_.clear();
    data_.clear();

    INT inactive=0;
    INT active=0;
    // loop over each time
    for(DataMapType::iterator element=measurements_.begin();element!=measurements_.end();++element)
      {
        bool activeMeasurements=false;
        // loop over each measurement
        for(size_t i=0;i<element->second.size();++i)
          {
            if (element->second[i].active)
              {
                element->second[i].index=data_.size();
                data_.push_back(element->second[i]);
                activeMeasurements=true;
                active++;
              }
            else
              inactive++;


            if (verbosity)
              std::cout << "value    " << left << fixed << setw(12) << setprecision(3) << element->second[i].value
                        << "type    " << left << fixed << setw(12)<< element->second[i].type
                        << "active    " << left << fixed << setw(4) << element->second[i].active
                        <<  "index    " << left << fixed << setw(4) << element->second[i].index
                        <<  "local i  " << left << fixed << setw(4) << element->second[i].localindex<< std::endl;
          }
        if (activeMeasurements)
          time_.push_back(element->first);

        if (verbosity)
          std::cout << "number of active measurements " << active << " and inactive " << inactive<< std::endl;
      }
  }



  void SetModelFileName(string filename)
  {
    modelFileName_=filename;
  }


  ~DycapFitDataClass()
  {}


private:
  DataMapType measurements_;
  vector<FLOAT> time_;
  vector<measurementType> data_;
  string modelFileName_;
  int outputCounter_;
  int numberofexperiments;
  int verbosity;
};

#endif
