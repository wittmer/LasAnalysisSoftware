///////////////////////////////////////////////////
// General use functions that act on std::vector //
///////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "LAS_vectorfloat_tools.h"

/////////////////////////////////////////////////////////////////////////////////////////
//! Print vector type and size to the output stream
std::ostream& operator<<(std::ostream& out, const std::vector<float> & svf)
{
  out << "std::vector<float> (size " << svf.size() << ")";
  return out;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Print vector type and size to the output stream
std::ostream& operator<<(std::ostream& out, const std::vector<double> & svd)
{
  out << "std::vector<double> (size " << svd.size() << ")";
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////
//! Fill the conetents of a vector into a histogram and plot it
TH1* vector_plot(const std::vector<double>& v, int nbins, const std::string& title, const std::string& name)
{
  //std::string hist_name = ((name == "") ? "vector_plot":name);

  // Create the limits for the x-axis
  double xmax=vector_max2(v);
  double xmin=vector_min2(v);
  double diff=xmax-xmin;
  xmax+= 0.1*diff;
  xmin-= 0.1*diff;
  std::cout << "xmin: " << xmin << "   xmax: " << xmax << std::endl;
  
  TH1 * h=new TH1D(name.c_str(),title.c_str(),nbins,xmin,xmax);
  
  //hfill(v,*h);
  for(std::vector<double>::size_type i=0; i<v.size(); i++){
    std::cout << "Filling " << v[i] << std::endl;
    h->Fill(v[i],1);
  }
  
  h->Draw();
  return h;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Draw two vectors of same size against each other
TGraph *vector_draw(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, const std::string & name)
{
  // Create new TGraph object;
  TGraph *gr = create_graph(v_x, v_y, name);
  if(gr == 0){
    std::cerr << "Error in vector_draw(v_x, v_y)" << std::endl;
    return 0;
  }

  gr->SetTitle(title.c_str());
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetTitle(xTitle.c_str());
  gr->GetYaxis()->SetTitle(yTitle.c_str());

  gr->Draw(options.c_str());

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Draw a single vector of floats in a TGraph object
TGraph *vector_draw(const std::vector<float> & data, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, const std::string & name)
{
  // Create new TGraph object;
  TGraph *gr = create_graph(data, name);

  gr->SetTitle(title.c_str());
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetTitle(xTitle.c_str());
  gr->GetYaxis()->SetTitle(yTitle.c_str());
  gr->Draw(options.c_str());

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Create a TGraph object from the vectors
TGraph *create_graph(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string & name)
{
  // Get the vector size
  std::vector<float>::size_type n_x = v_x.size();
  std::vector<float>::size_type n_y = v_y.size();
  if(n_x != n_y){
    std::cerr << "Error in create_graph(v_x, v_y): The vecors do not have the same size! v_x.size(): " << v_x.size() << "  v_y.size(): " << v_y.size() << std::endl;
    return 0;
  }

  // Create new TGraph object;
  TGraph *gr = new TGraph(n_x);
  gr->SetName(name.c_str());

  // Fill in the data
  for(std::vector<float>::size_type i=0; i<n_x; i++)
    gr->SetPoint(i, v_x[i] , v_y[i]);
  gr->Sort();

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Create a TGraphErrors object from the vectors
TGraphErrors* create_graph_errors(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string & name)
{
  // Get the vector size
  std::vector<float>::size_type n_x = v_x.size();
  std::vector<float>::size_type n_y = v_y.size();
  if(n_x != n_y){
    std::cerr << "Error in create_graph_errors(v_x, v_y): The vecors do not have the same size! v_x.size(): " << v_x.size() << "  v_y.size(): " << v_y.size() << std::endl;
    return 0;
  }

  // Create new TGraph object;
  TGraphErrors *gr = new TGraphErrors(n_x);
  gr->SetName(name.c_str());

  // Fill in the data
  for(std::vector<float>::size_type i=0; i<n_x; i++){
    gr->SetPoint(i, v_x[i] , v_y[i]);
    gr->SetPointError(i, 0, 2);
  }
  gr->Sort();

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Create a TGraph object and fill it with the data values. The graph will not be drawn
TGraph *create_graph(const std::vector<float> & data, const std::string & name)
{
  // Get the vector size
  std::vector<float>::size_type n = data.size();

  // Create new TGraph object;
  TGraph *gr = new TGraph(n);
  gr->SetName(name.c_str());

  // Fill in the data
  for(std::vector<float>::size_type i=0; i<n; i++)
    gr->SetPoint(i, i+1 , data[i]);
  gr->Sort();

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Create a TGraphErrors object and fill it with the data values. The graph will not be drawn
TGraphErrors *create_graph_errors(const std::vector<float> & data, const std::string & name, double error)
{
  // Get the vector size
  std::vector<float>::size_type n = data.size();

  // Create new TGraph object;
  TGraphErrors *gr = new TGraphErrors(n);
  gr->SetName(name.c_str());

  // Fill in the data
  for(std::vector<float>::size_type i=0; i<n; i++){
    gr->SetPoint(i, i+1 , data[i]);
    gr->SetPointError(i, 0, error);
  }
  gr->Sort();

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Compute mean and rms of the vector entries
void average_positions(const std::vector<double> data, double& mean, double& rms)
{
  if(data.empty()){
    std::cerr << "Error in average_positions: The data is empty" << std::endl;
    return;
  }
  double sum=0;
  double sumsquares = 0;
  for(std::vector<double>::size_type i = 0 ; i < data.size(); i++){
    sum += data[i];
    sumsquares += data[i]*data[i];
  }
  mean = sum / data.size();
  rms = sqrt(sumsquares/data.size() - mean*mean);
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Add all entries from data to sum
void vector_add(std::vector<float>& sum, const std::vector<float>& data)
{
  std::vector<float>::size_type size = sum.size();
  if(size != data.size()){
    std::cerr << "Error in vector_add, the two vectors do not have the same size. data.size(): " << data.size() << "   sum.size(): " << sum.size() << std::endl;
    return;
  }

  for(std::vector<float>::size_type i = 0; i < size; i++) sum[i] += data[i];
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Divide al entries in the data vector by norm
void vector_norm(std::vector<float>& data, float norm)
{
  if(norm == 0){
    std::cerr << "Error in vector_norm, the normalization factor is zero!" << std::endl;
    return;
  }

  for(std::vector<float>::size_type i = 0; i < data.size(); i++) data[i] /= norm;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the maximum of the vector entries (range from low to high)
float vector_max(const std::vector<float>& v, vfloat_idx low, vfloat_idx high)
{
  // First do some basic checks and range adjustments
  vfloat_idx v_lo = (low < 0) ? 0          :  low;
  vfloat_idx v_hi = (high< 0 || high >= v.size()) ? v.size()-1 : high;
  if( v.empty() || v_hi <= v_lo )return 0;

  return *(std::max_element(v.begin() + v_lo, v.begin() + v_hi + 1)); 
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the maximum of the vector entries (range from low to high)
float vector_max_old(const std::vector<float>& v, int low, int high)
{
  // First do some basic checks and range adjustments
  vfloat_idx v_lo = (low < 0) ? 0          :  low;
  vfloat_idx v_hi = (high< 0) ? v.size()-1 : high;
  if( v.empty() || v_lo >= v.size() )return 0;
  if(v_hi < v_lo || v_hi >= v.size()) v_hi = (v.size() - 1);
 
  float max = v[v_lo];

  for( vfloat_idx idx = v_lo; idx <= v_hi; idx++){
    max = std::max( max , v[idx]);
  }
  return max;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the index with the highest entry (range from low to high)
vfloat_idx vector_max_idx(const std::vector<float>& v, int low, int high)
{
  // First do some basic checks and range adjustments
  vfloat_idx v_lo = (low < 0) ? 0          :  low;
  vfloat_idx v_hi = (high< 0) ? v.size()-1 : high;
  if( v.empty() || v_lo >= v.size() )return 0;
  if(v_hi < v_lo || v_hi >= v.size()) v_hi = (v.size() - 1);
 
  int max_idx= v_lo;

  for( vfloat_idx idx = v_lo; idx <= v_hi; idx++){
    if( v[max_idx] < v[idx]) max_idx = idx ;
  }
  return max_idx;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the maximum of the vector entries (this version has slightly better performace but no ranges)
float vector_max2(const std::vector<float>& v)
{
  if(v.empty())return 0;
  std::vector<float>::const_iterator it = v.begin();
  float max = *it;
  while(it++ != v.end()){
    max = std::max( max , *it);
  }
  return max;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the maximum of the vector entries (this version has slightly better performace but no ranges)
double vector_max2(const std::vector<double>& v)
{
  if(v.empty())return 0;
  std::vector<double>::const_iterator it = v.begin();
  double max = *it;
  while(it++ != v.end()){
    max = std::max( max , *it);
  }
  return max;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the minmum of the vector entries (this version has slightly better performace but no ranges)
float vector_min2(const std::vector<float>& v)
{
  if(v.empty())return 0;
  std::vector<float>::const_iterator it = v.begin();
  float min = *it;
  while(it++ != v.end()){
    min = std::min( min , *it);
  }
  return min;
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Return the minmum of the vector entries (this version has slightly better performace but no ranges)
float vector_min2(const std::vector<double>& v)
{
  if(v.empty())return 0;
  std::vector<double>::const_iterator it = v.begin();
  double min = *it;
  while(it++ != v.end()){
    min = std::min( min , *it);
  }
  return min;
}
