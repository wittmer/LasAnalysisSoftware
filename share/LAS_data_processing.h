#ifndef __LAS_DATA_BROCESSING_H__
#define __LAS_DATA_BROCESSING_H__

void main_analysis(const std::string& filename, const std::string& output_filename, Long64_t nentries = -1, bool max_signal = true, bool noise_filter = true, bool block_slicer = true, bool labels = true, bool profiles = true, bool positions = true, bool tec_reco = true);

void history_plots(const std::string& listfile, const std::string& outfile, Bool_t suppress_graph = kTRUE);

void compress_AT_plots(const std::string& parameters_file, Avec::size_type n = 100, const std::string& directory = "Details_AT/data/", LAS::beam_group beam_group = LAS::AT);


#endif
