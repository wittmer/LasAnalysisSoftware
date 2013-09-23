#ifndef __LAS_VECTORFLOAT_TOOLS__
#define __LAS_VECTORFLOAT_TOOLS__

#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// Typedefs to make typing easier
typedef std::vector<float>::size_type vfloat_idx;
typedef std::vector<double>::size_type vdouble_idx;

//////////////////////////////
TGraph *vector_draw(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, const std::string & name);
TGraph *vector_draw(const std::vector<float> & data, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, const std::string & name);
TGraph* create_graph(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string & name);
TGraph* create_graph(const std::vector<float> & data, const std::string & name);
TGraphErrors *create_graph_errors(const std::vector<float> & data, const std::string & name, double error);
TGraphErrors* create_graph_errors(const std::vector<float> & v_x, const std::vector<float> & v_y, const std::string & name);
TH1* vector_plot(const std::vector<double>& v, int nbins, const std::string& title, const std::string& name);

void average_positions(const std::vector<double> data, double& mean, double& rms);
void vector_add(std::vector<float>& sum, const std::vector<float>& data);
void vector_norm(std::vector<float>& data, float norm);

float vector_max(const std::vector<float>& v, vfloat_idx low=0, vfloat_idx high=-1);
vfloat_idx vector_max_idx(const std::vector<float>& v, int low=0, int high=-1);

float vector_max2(const std::vector<float>& v);
double vector_max2(const std::vector<double>& v);
float vector_min2(const std::vector<float>& v);
float vector_min2(const std::vector<double>& v);

#endif
