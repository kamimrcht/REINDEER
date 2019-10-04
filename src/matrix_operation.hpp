#ifndef MAT
#define MAT





using namespace std;


	


void write_color_matrix(const string& output_file, vector<vector<uint8_t>>& color_matrix){
	auto out=new zstr::ofstream(output_file);
	uint64_t color_number(color_matrix.size());
	uint64_t line_size(color_matrix[0].size());
	uint i(0);

	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t));
	out->write(reinterpret_cast<char*>(&line_size),sizeof(uint64_t));
	//~ auto out=new zstr::ofstream(output_file);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		auto point =&(color_matrix[vec][0]);
		out->write((char*)point,(line_size));
	}
	delete out;
}
void write_color_matrix_counts(const string& output_file, vector<vector<uint16_t>>& color_matrix){
	auto out=new zstr::ofstream(output_file);
	uint64_t color_number(color_matrix.size());
	uint64_t line_size(color_matrix[0].size());
	uint i(0);

	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t));
	out->write(reinterpret_cast<char*>(&line_size),sizeof(uint64_t));
	//~ auto out=new zstr::ofstream(output_file);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		auto point =&(color_matrix[vec][0]);
		out->write((char*)point,(line_size));
	}
	delete out;
}

void write_color_matrix_reads(const string& output_file, vector<vector<uint32_t>>& color_matrix){
	auto out=new zstr::ofstream(output_file);
	uint64_t color_number(color_matrix.size());
	uint64_t line_size(color_matrix[0].size());
	uint i(0);

	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t));
	out->write(reinterpret_cast<char*>(&line_size),sizeof(uint64_t));
	//~ auto out=new zstr::ofstream(output_file);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		auto point =&(color_matrix[vec][0]);
		out->write((char*)point,(line_size));
	}
	delete out;
}



vector<vector<uint8_t>> load_written_matrix(const string& input_file){
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	uint64_t color_number;
	uint64_t line_size;
	auto in=new zstr::ifstream(input_file);
	in-> read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	in-> read(reinterpret_cast<char *>(&line_size), sizeof(uint64_t));
	vector<uint8_t> colorV(line_size,0);
	vector<vector<uint8_t>> color_matrix(color_number,colorV);
	//~ auto in=new zstr::ifstream(input);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		in->read( (char*)((color_matrix[vec].data())) , line_size);
	}
	uint i(0);
	delete(in);
	return color_matrix;
}
vector<vector<uint16_t>> load_written_matrix_counts(const string& input_file){
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	uint64_t color_number;
	uint64_t line_size;
	auto in=new zstr::ifstream(input_file);
	in-> read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	in-> read(reinterpret_cast<char *>(&line_size), sizeof(uint64_t));
	vector<uint16_t> colorV(line_size,0);
	vector<vector<uint16_t>> color_matrix(color_number,colorV);
	//~ auto in=new zstr::ifstream(input);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		in->read( (char*)((color_matrix[vec].data())) , line_size);
	}
	uint i(0);
	delete(in);
	return color_matrix;
}

vector<vector<uint32_t>> load_written_matrix_reads(const string& input_file){
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	uint64_t color_number;
	uint64_t line_size;
	auto in=new zstr::ifstream(input_file);
	in-> read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	in-> read(reinterpret_cast<char *>(&line_size), sizeof(uint64_t));
	vector<uint32_t> colorV(line_size,0);
	vector<vector<uint32_t>> color_matrix(color_number,colorV);
	//~ auto in=new zstr::ifstream(input);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		in->read( (char*)((color_matrix[vec].data())) , line_size);
	}
	uint i(0);
	delete(in);
	return color_matrix;
}



void dump_map(const string& output_file_vector, const string& output_file_position, const string& output_file_size,  vector<vector<uint64_t>>& unitigs_to_nodes)
{
	auto out_vector=new ofstream(output_file_vector);
	auto out_position=new ofstream(output_file_position);
	auto out_size=new ofstream(output_file_size);
	for (auto && unitig: unitigs_to_nodes){
		uint32_t vec_size(unitig.size());
		out_size->write(reinterpret_cast<char*>(&vec_size),sizeof(uint32_t));
		uint64_t file_pos(out_vector->tellp()/8);
		out_position->write(reinterpret_cast<char*>(&file_pos),sizeof(uint64_t));
		if (vec_size > 0)
		{
			auto point =&(unitig[0]);
			out_vector->write((char*)point,(vec_size*sizeof(uint64_t)));
		}
	}
	delete out_vector;
	delete out_position;
	delete out_size;
}

void load_unitig(uint32_t unitigID, const string& file_vector, const string& file_position, const string& file_size, vector<uint64_t>& reads)
{
	uint32_t vec_size(0);
	auto in_size=new ifstream(file_size);
	in_size->seekg(unitigID*sizeof(uint32_t), in_size->beg);
	in_size-> read(reinterpret_cast<char *>(&vec_size), sizeof(uint32_t));
	if (vec_size > 0){
		uint64_t vec_position;
		auto in_position=new ifstream(file_position);
		in_position->seekg(unitigID*sizeof(uint64_t), in_position->beg);
		in_position-> read(reinterpret_cast<char *>(&vec_position), sizeof(uint64_t));
		auto in_vector=new ifstream(file_vector);
		in_vector->seekg(vec_position*sizeof(uint64_t), in_vector->beg);
		reads.resize(vec_size,0);
		in_vector-> read(reinterpret_cast<char *>(reads.data()), vec_size*sizeof(uint64_t));
		delete(in_position);
		delete(in_vector);
	}
	delete(in_size);
}



#endif
