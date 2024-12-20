#ifndef QUERYFUNCTION_H
#define QUERYFUNCTION_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "DSTGS.h"
using namespace std;

const int query_data_pairs = 100000;

struct DSTGS_VAR { 
	time_type startTime;
	uint32_t granularityLength;
	uint32_t gl;
	uint32_t width;
	uint32_t depth;
	uint32_t fingerprintLen;
	uint32_t row_addrs;
	uint32_t col_addrs;
};

struct QueryPairData {
	uint64_t source;
	uint64_t destination;
	time_type start_time;
	time_type end_time;
};

/******************* static variables *************************/
static DSTGS* dstgs;
/******************* static variables *************************/
/***************** function declaration ***********************/
// basic functions
int isFolderExist(char* folder);
int createDirectory(char* sPathName);
uint64_t count_lines(string file);
int readRandomFileToDataArray(string file, QueryPairData dataArray[]);
// insert function
int baselineInsert(DSTGS_VAR var, string filename);
// para query functions
int64_t edgeFrequenceBaseline(DSTGS& DSTGS, int64_t s, int64_t d, int64_t start, int64_t end);
int64_t nodeFrequenceBaseline(DSTGS& DSTGS, int64_t v, int type, int64_t start, int64_t end);
int edgeFrequenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write);
int edgeFrequenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write);
int edgeExistenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write, int flag);
int edgeExistenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag);
int nodeFrequenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write, int flag, int line);
int nodeFrequenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line);
// seq query functions
int edgeFrequenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write);
int edgeExistenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag);
int nodeFrequenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line);
// query functions that the main function called
void edgeFrequenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write);
void edgeExistenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag);
void nodeFrequenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line);
/***************** function declaration ***********************/

int isFolderExist(char* folder) {
	int ret = 0;
	ret = access(folder, R_OK | W_OK);
	if (ret == 0)
		ret = 1;
	else
		ret = 0;
	return ret;
}
int createDirectory(char* sPathName) {
	char DirName[256];
	strcpy(DirName, sPathName);
	int i, len = strlen(DirName);
	if (DirName[len - 1] != '/')
		strcat(DirName, "/");
	len = strlen(DirName);
	for (i = 1; i < len; i++) {
		if (DirName[i] == '/') {
			DirName[i] = 0;
			int a = access(DirName, F_OK);
			if (a == -1) {
				mkdir(DirName, 0755);
			}
			DirName[i] = '/';
		}
	}
	return 0;
}
uint64_t count_lines(string file) {
	ifstream readFile;
	uint64_t n = 0;
	char line[512];
	string temp;
	readFile.open(file, ios::in);	// ios::in means that open file with readonly 
	if(readFile.fail()) { 			// open file error, return 0
		cout << "error in opening file" << endl;
		return 0;
	}
	else { 							// the file exists
		while(getline(readFile,temp))
			n++;
	}
	readFile.close();
	return n;
}
int readRandomFileToDataArray(string file, QueryPairData dataArray[]) {
	ifstream randomFile;
	randomFile.open(file);
	if (!randomFile.is_open()) {
		cout << "Error in open file, Path = " << file << endl;
		return -1;
	}
	int datanum = 0;
	int64_t startPoint, endPoint, timeStart, timeEnd;
	while (!randomFile.eof()) {
		randomFile >> startPoint >> endPoint >> timeStart >> timeEnd;
		if(randomFile.fail())
			break;
		dataArray[datanum].source = startPoint;
		dataArray[datanum].destination = endPoint;
		dataArray[datanum].start_time = timeStart;
		dataArray[datanum].end_time = timeEnd;
		datanum++;
		// if(datanum > query_data_pairs) {
		// 	cout << "the input data is more than the range of the array" << endl;
		// 	break;
		// }
	}
	randomFile.close();
	return datanum;
}
#if defined(DEBUG) || defined(TINSTIME) || defined(BINSTIME) || defined(HINT)
void progress_bar(int n) {
	int i = 0;
	char bar[102];
	const char *lable = "|/-\\";
	bar[0] = 0;
	while (i < n) {
		bar[i] = '#';
		i++;
		bar[i] = 0;
	}
	printf("\r[%-100s][%d%%][%c]", bar, i, lable[i%4]);
	fflush(stdout);
	return;
}
#endif
//边插入
int baselineInsert(DSTGS_VAR var, string filename) {
    //1.打开数据集文件
	ifstream ifs;
	ifs.open(filename);
	if (!ifs.is_open()) {
		cout << "Error in open file, Path = " << filename << endl;
		return -1;
	}
	int64_t s, d;
	weight_type w;
	time_type t;
#if defined(DEBUG) || defined(TINSTIME) || defined(BINSTIME) || defined(HINT)
	cout << "Inserting..." << endl;
	timeval matrix_s, matrix_e;
	gettimeofday( &matrix_s, NULL);
#endif
	//2.创建DSTGS对象（调用DSFGS的构造函数）
	dstgs = new DSTGS (var.gl, var.width, var.depth, var.fingerprintLen, var.row_addrs, var.col_addrs);//var.gl=1
	cout << "创建DSTGS成功" <<endl;

#if defined(DEBUG) || defined(TINSTIME)
	gettimeofday( &matrix_e, NULL);
	double matrix_time = (matrix_e.tv_sec - matrix_s.tv_sec) +  (matrix_e.tv_usec - matrix_s.tv_usec) / 1000000.0;
	cout << "Matrix Time = " << matrix_time << " s" << endl;
#endif
#if defined(DEBUG) || defined(TINSTIME) || defined(BINSTIME) || defined(HINT)
	timeval ins_start, ins_end;	
	gettimeofday( &ins_start, NULL);
#endif
#if defined(DEBUG) || defined(BINSTIME)
	timeval bpoint_start, bpoint_end;
	gettimeofday( &bpoint_start, NULL);
#endif
#if defined(DEBUG) || defined(HINT)
		double total = count_lines(filename);
		if(total == 0)
			cout << "ERROR--QueryFunction--178" <<endl;
#endif
    //3.开始读取文件中的每条流边
	cout << "开始读取流边" <<endl;
	int datanum = 0;
	while (!ifs.eof()) {
		ifs >> s >> d >> w >> t;
		if(ifs.fail())
			break;
        //计算当前流边的时间段tt
		uint32_t tt = ceil((double)(t - var.startTime) / (double)var.granularityLength);
        //调用dstgs的insert方法，完成当前边的插入，传入的参数为源节点s，目的节点d，边权重w和流边所属时间段tt
		string sv=to_string(s);
		string dv=to_string(d);
		dstgs->insert(sv, dv, w, tt);
		datanum++;

#if defined(DEBUG) || defined(BINSTIME)
		if (datanum % 10000000 == 0) {
			gettimeofday( &bpoint_end, NULL);
			double bpoint_time = (bpoint_end.tv_sec - bpoint_start.tv_sec) +  (bpoint_end.tv_usec - bpoint_start.tv_usec) / 1000000.0;
			cout << datanum << ", Break Point Time = " << bpoint_time << " s" << endl;
			gettimeofday( &bpoint_start, NULL);
		}
#endif
#if defined(DEBUG) || defined(HINT)
		if (datanum % 100000 == 0) {
			int n = (int) ((double) datanum / total * 100);
			progress_bar(n);
		}
		if (datanum == total) {
			progress_bar(100);
		}
#endif
	}
    //流边插入结束
#if defined(DEBUG) || defined(HINT)
	cout << endl;
#endif
#if defined(DEBUG) || defined(BINSTIME)
	gettimeofday( &bpoint_end, NULL);
	double bpoint_time = (bpoint_end.tv_sec - bpoint_start.tv_sec) +  (bpoint_end.tv_usec - bpoint_start.tv_usec) / 1000000.0;
	cout << datanum << ", Break Point Time = " << bpoint_time << " s" << endl;
#endif
#if defined(DEBUG) || defined(TINSTIME) || defined(BINSTIME) || defined(HINT)
	gettimeofday( &ins_end, NULL);
	double ins_time = (ins_end.tv_sec - ins_start.tv_sec) +  (ins_end.tv_usec - ins_start.tv_usec) / 1000000.0;
	cout << "Insertion Time = " << ins_time << " s" << endl;
	cout << "Insertion Finished!" << endl;
	cout << "Datanum = " << datanum << endl;
#endif


#if defined(DEBUG) || defined(TINSTIME) || defined(BINSTIME) || defined(HINT)
	cout << "Flashing..." << endl;
	timeval flash_s, flash_e;
	gettimeofday( &flash_s, NULL);
#endif
	//开启flash
	dstgs->flash();
#if defined(DEBUG) || defined(TINSTIME)|| defined(BINSTIME) || defined(HINT)
	gettimeofday( &flash_e, NULL);
	double flash_time = (flash_e.tv_sec - flash_s.tv_sec) +  (flash_e.tv_usec - flash_s.tv_usec) / 1000000.0;
	cout << "Flash Time = " << flash_time << " s" << endl;
#endif

	cout << "************************" << endl;
	dstgs->bucketCounting();
	cout << "************************" << endl << endl;
	//计算并输出内存大小
	dstgs->updateMemorySizes();
    if (dstgs!= nullptr)
    {
        size_t totalMemorySize = dstgs->getTotalMemorySize();
        cout << "插入操作结束后使用的内存大小为: " << totalMemorySize << " 字节" << endl;
    }
    else
    {
        cout << "插入操作可能未正确执行，无法获取内存使用大小" << endl;
    }
	ifs.close();
	return 0;
}
//并行查询（具体流程：para调用single，single调用baseline）
int64_t edgeFrequenceBaseline(DSTGS& DSTGS, int64_t s, int64_t d, int64_t start, int64_t end) {
	string sv=to_string(s);
	string dv=to_string(d);
	// int64_t result = 0;
    // //目标查询范围分解
    // vector<interval> subranges=decomposeTargetRange(start,end);
    // vector<interval> subintervals=furtherDecompose(subranges);
	// cout << subintervals.size() << endl;
    // //对每一个可查询子区间执行边查询
    // for (const auto& sub : subintervals) {
    //     result += DSTGS.edgeQuery(sv,dv,sub.start,sub.end);
    // }
	// return result;
	int64_t result = 0;
	result +=DSTGS.edgeQuery(sv,dv,start,end);
	return result;
}
int64_t nodeFrequenceBaseline(DSTGS& DSTGS, int64_t v, int type, int64_t start, int64_t end){
	string vv=to_string(v);
    int64_t result = 0;
    //目标查询范围分解
	//vector<interval> subintervals=decomposeTargetRange(start,end);
    //对每一个可查询子区间执行点查询
    // for (const auto& sub : subintervals) {
    //     result += DSTGS.nodeQuery(vv,type,sub.start,sub.end);
    // }
	result +=DSTGS.nodeQuery(vv,type,start,end);
	return result;
}
int edgeFrequenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）
    string input_file_prefix = dataset_name + "_random_edge_frequence_";
	string input_file_suffix = "_sorted.txt";
	string output_file_prefix = dataset_name + "_edge_frequence_baseline_";
	string output_file_suffix = "_res.txt";
	string time_file_suffix = "_time.txt";

	//2.读取测试集（单个）中的每一个查询项，记入dataArray中，共datanum项
	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
	int datanum = readRandomFileToDataArray(input_dir + input_file_prefix + to_string(num) + input_file_suffix, dataArray);

#if defined(DEBUG) || defined(HINT)
	cout << "****************** timeslot = " << num << " ******************" << endl;
#endif

    //3.创建输出文件
	ofstream resultFile, timeFile;
	if (write) {
		char dir_path[FILENAME_MAX];
		strcpy(dir_path, output_dir.c_str());
		if (createDirectory(dir_path) != 0) {
			cout << "CreateDirectory Error, Path = " << dir_path << endl;
			return -1;
		}
		resultFile.open(output_dir + output_file_prefix + to_string(num) + output_file_suffix);
		if (!resultFile.is_open()) {
			cout << "Error in open file, Path = " << (output_dir + output_file_prefix + to_string(num) + output_file_suffix) << endl;
			return -1;
		}
		timeFile.open(output_dir + output_file_prefix + to_string(num) + time_file_suffix);
		if (!timeFile.is_open()) {
			cout << "Error in open file, Path = " << (output_dir + output_file_prefix + to_string(num) + time_file_suffix) << endl;
			return -1;
		}
	}

    //4.开始对每一个查询项执行边权重查询
	double sumTime = 0, sumTime_perquery = 0;
	timeval tp1, tp2;
	for (int m = 0; m < query_times; m++) {
		sumTime_perquery = 0;
		for (int n = 0; n < datanum; n++) {
			gettimeofday( &tp1, NULL);
            //每一个查询项都调用edgeFrequenceBaseline方法
			int64_t res = edgeFrequenceBaseline(*dstgs, dataArray[n].source, dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
			gettimeofday( &tp2, NULL);
			double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
			sumTime_perquery += delta_t;
			if (write && m == 0) {
				if(n == (datanum - 1)) {
					resultFile << res;
					timeFile  << delta_t;
					break;
				}
				else {
					resultFile << res << endl;
					timeFile  << delta_t << endl;
				}
			}
		}
		sumTime += (sumTime_perquery / (double)datanum);
	}

	if (write) {
		resultFile.flush();
		timeFile.flush();
		resultFile.close();
		timeFile.close();
	}
	delete[] dataArray;

#if defined(DEBUG) || defined(HINT)
	double mseconds = (double)(sumTime / (double)query_times) / 1000;
	printf("Timeslot = %d, Query Times = %d, Query Avg Time = %lf ms\n\n", num, query_times, mseconds);
#endif

	return 0;
}
int edgeFrequenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write){
    vector<thread*> childs;
	thread* child = NULL;
	for(int i = 0; i < num.size(); i++) {
		child = new thread(edgeFrequenceBaselineTest_single, input_dir, output_dir, dataset_name, num[i], query_times, write);
		childs.push_back(child);

	}
	for(int i = 0; i < childs.size(); i++) {
		if(childs[i] != NULL)
			childs[i]->join();
	}
	return 0;
}
int edgeExistenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write, int flag){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）,flag仅在确定输入和输出文件上起作用，其他执行过程完全相同
    string input_file_prefix = "";
	string input_file_suffix = "";
	string output_file_prefix = "";
	string output_file_suffix = "";
	string time_file_suffix = "_time.txt";
	switch (flag)
	{
	case 1:
		input_file_prefix = "_random_edge_existence_";
		input_file_suffix = "_sorted.txt";
		output_file_prefix = "_edge_existence_baseline_";
		output_file_suffix = "_res.txt";
		break;
	case 2:
		input_file_prefix = "_bool_";
		input_file_suffix = ".txt";
		output_file_prefix = "_bool_baseline_";
		output_file_suffix = "_res.txt";
		break;
	default:
		break;
	}

	//2.读取测试集（单个）中的每一个查询项，记入dataArray中，共datanum项
	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
	int datanum = readRandomFileToDataArray(input_dir + dataset_name + input_file_prefix + to_string(num) + input_file_suffix, dataArray);

#if defined(DEBUG) || defined(HINT)
	cout << "****************** timeslot = " << num << " ******************" << endl;
#endif

    //3.创建输出文件
	ofstream resultFile, timeFile;
	if (write) {
		char dir_path[FILENAME_MAX];
		strcpy(dir_path, output_dir.c_str());
		if (createDirectory(dir_path) != 0) {
			cout << "createDirectory error" << endl;
			return -1;
		}
		resultFile.open(output_dir + dataset_name + output_file_prefix + to_string(num) + output_file_suffix);
		if (!resultFile.is_open()) {
			cout << "error in open file " << (output_dir + dataset_name + output_file_prefix + to_string(num) + output_file_suffix) << endl;
			return -1;
		}
		timeFile.open(output_dir + dataset_name + output_file_prefix + to_string(num) + time_file_suffix);
		if (!timeFile.is_open()) {
			cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + to_string(num) + time_file_suffix) << endl;
			return -1;
		}
	}

    //4.开始对每一个查询项执行边权重查询
	double sumTime = 0, sumTime_perquery = 0;
	int ones = 0;
	timeval tp1, tp2;
	for (int m = 0; m < query_times; m++) {
		sumTime_perquery = 0;
		for (int n = 0; n < datanum; n++) {
			gettimeofday( &tp1, NULL);
            //每个查询项都调用edgeFrequenceBaseline方法
			int64_t res = edgeFrequenceBaseline(*dstgs, dataArray[n].source, dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
			gettimeofday( &tp2, NULL);
			double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
			sumTime_perquery += delta_t;
			if (write && m == 0) {
				if (res > 0)   
					ones++;
				if(n == (datanum - 1)) {
					resultFile << ((res > 0) ? 1 : 0);
					timeFile  << delta_t;
					break;
				}
				else {
					resultFile << ((res > 0) ? 1 : 0) << endl;
					timeFile  << delta_t << endl;
				}
			}
		}
		sumTime += (sumTime_perquery / (double)datanum);
	}
	
	if (write) {
		resultFile.flush();
		timeFile.flush();
		resultFile.close();
		timeFile.close();
	}
	delete[] dataArray;
#if defined(DEBUG) || defined(HINT)
	double mseconds = (double)(sumTime / (double)query_times) / 1000; 
	printf("Timeslot = %d, Query Times = %d, Query Avg Time = %lf ms\n\n", num, query_times, mseconds);
#endif
	return 0;
}
int edgeExistenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag){
    vector<thread*> childs;
	thread* child = NULL;
	for(int i = 0; i < num.size(); i++) {
		child = new thread(edgeExistenceBaselineTest_single, input_dir, output_dir, dataset_name, num[i], query_times, write, flag);
		childs.push_back(child);

	}
	for(int i = 0; i < childs.size(); i++) {
		if(childs[i] != NULL)
			childs[i]->join();
	}
	return 0;
}
int nodeFrequenceBaselineTest_single(string input_dir, string output_dir, string dataset_name, int num, int query_times, bool write, int flag, int line){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）,flag仅在确定输入和输出文件上起作用（在后面的执行过程中，不同文件中的第二个元素发挥type的作用）
    string input_file_prefix = "";
	string input_file_suffix = "";
	string output_file_prefix = "";
	string output_file_suffix = "";
	string time_file_suffix = "_time.txt";
	switch (flag) {
		case 1:
			input_file_prefix = "_random_node_frequence_in_";
			input_file_suffix = "_sorted.txt";
			output_file_prefix = "_node_frequence_in_";
			output_file_suffix = "_res.txt";
			break;
		case 2:
			input_file_prefix = "_random_node_frequence_out_";
			input_file_suffix = "_sorted.txt";
			output_file_prefix = "_node_frequence_out_";
			output_file_suffix = "_res.txt";
			break;
		default:
			break;
	}

	//2.读取测试集（单个）中的每一个查询项，记入dataArray中，共datanum项
	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
	int datanum = readRandomFileToDataArray(input_dir + dataset_name + input_file_prefix + to_string(num) + input_file_suffix, dataArray);

#if defined(DEBUG) || defined(HINT)
	cout << "****************** timeslot = " << num << " ******************" << endl;
#endif

    //3.创建输出文件
	ofstream resultFile, timeFile;
	if (write) {
		char dir_path[FILENAME_MAX];
		strcpy(dir_path, output_dir.c_str());
		if (createDirectory(dir_path) != 0) {
			cout << "createDirectory error" << endl;
			return -1;
		}

		if(line == 0) {
			resultFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + output_file_suffix);
			if (!resultFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + output_file_suffix) << endl;
				return -1;
			}
			timeFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + time_file_suffix);
			if (!timeFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + time_file_suffix) << endl;
				return -1;
			}
		}
		else {	//	append
			resultFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + output_file_suffix, ios::app);
			cout << "append " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + output_file_suffix) << endl;
			if (!resultFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + output_file_suffix) << endl;
				return -1;
			}
			timeFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + time_file_suffix, ios::app);
			cout << "append " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + time_file_suffix) << endl;
			if (!timeFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num) + time_file_suffix) << endl;
				return -1;
			}
		}
	}

    //4.开始对每一个查询项执行点权重查询
	double sumTime = 0, sumTime_perquery = 0;	
	timeval tp1, tp2;
	for (int m = 0; m < query_times; m++) {
		sumTime_perquery = 0;
		for (int n = 0; n < datanum; n++) {
			if(m == 0 && n < line) 
				continue;		
			gettimeofday( &tp1, NULL);
            //每个查询项调用nodeFrequenceBaseline方法
			int64_t res = nodeFrequenceBaseline(*dstgs, dataArray[n].source, (int)dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
			gettimeofday( &tp2, NULL);
			double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
			sumTime_perquery += delta_t;
			
			if (write && m == 0) {
				if(n == (datanum - 1)) {
					resultFile << res;
					timeFile  << delta_t;
					break;
				}
				else {
					resultFile << res << endl;
					timeFile  << delta_t << endl;
				}
			}
		}
		sumTime += (sumTime_perquery / (double)datanum);
	}

	if (write) {
		resultFile.flush();
		timeFile.flush();
		resultFile.close();
		timeFile.close();
	}
	delete[] dataArray;
#if defined(DEBUG) || defined(HINT)
	double mseconds = (double)(sumTime / (double)query_times) / 1000; 
	printf("Timeslot = %d, Query Times = %d, Query Avg Time = %lf ms\n\n", num, query_times, mseconds);
#endif
	return 0;
}
int nodeFrequenceBaselineTest_para(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line){
    vector<thread*> childs;
	thread* child = NULL;
	for(int i = 0; i < num.size(); i++) {
		child = new thread(nodeFrequenceBaselineTest_single, input_dir, output_dir, dataset_name, num[i], query_times, write, flag, line);
		childs.push_back(child);

	}
	for(int i = 0; i < childs.size(); i++) {
		if(childs[i] != NULL)
			childs[i]->join();
	}
	return 0;
}
//顺序查询（具体流程，seq调用baseline）
int edgeFrequenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）
    string input_file_prefix = dataset_name + "_random_edge_frequence_";
	string input_file_suffix = "_sorted.txt";
	string output_file_prefix = dataset_name + "_edge_frequence_baseline_";
	string output_file_suffix = "_res.txt";
	string time_file_suffix = "_time.txt";

	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
    //2.一个一个测试文件按顺序遍历、读取数据项、执行边权重查询
	for (int i = 0; i < num.size(); i++) {
        //2.1读取当前测试文件内容到dataArray
		int datanum = readRandomFileToDataArray(input_dir + input_file_prefix + to_string(num[i]) + input_file_suffix, dataArray);

#if defined(DEBUG) || defined(HINT)
		cout << "****************** timeslot = " << num[i] << " ******************" << endl;
#endif
        //2.2创建当前测试文件对应的输出文件
		ofstream resultFile, timeFile;
		if (write) {
			char dir_path[FILENAME_MAX];
			strcpy(dir_path, output_dir.c_str());
			if (createDirectory(dir_path) != 0) {
				cout << "CreateDirectory Error, Path = " << dir_path << endl;
				return -1;
			}
			resultFile.open(output_dir + output_file_prefix + to_string(num[i]) + output_file_suffix);
			if (!resultFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + output_file_prefix + to_string(num[i]) + output_file_suffix) << endl;
				return -1;
			}
			timeFile.open(output_dir + output_file_prefix + to_string(num[i]) + time_file_suffix);
			if (!timeFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + output_file_prefix + to_string(num[i]) + output_file_suffix) << endl;
				return -1;
			}
		}
        //2.3开始对当前dataArray中的每一个查询项执行边权重查询
		double sumTime = 0, sumTime_perquery = 0;
		timeval tp1, tp2;
		for (int m = 0; m < query_times; m++) {
			sumTime_perquery = 0;
			for (int n = 0; n < datanum; n++) {
				gettimeofday( &tp1, NULL);
				//对文件当中的每条查询调用edgeFrequenceBaseline方法，传入源节点、目的节点、起始时间段、结束时间段
				int64_t res = edgeFrequenceBaseline(*dstgs, dataArray[n].source, dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
				gettimeofday( &tp2, NULL);
				double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
				sumTime_perquery += delta_t;
				if (write && (m == 0)) {
					if(n == (datanum - 1)) {
						resultFile << res;
						timeFile  << delta_t;
						break;
					}
					else {
						resultFile << res << endl;
						timeFile  << delta_t << endl;
					}
				}
			}
			sumTime += (sumTime_perquery / (double)datanum);
		}
		if (write) {
			resultFile.flush();
			timeFile.flush();
			resultFile.close();
			timeFile.close();
		}
#if defined(DEBUG) || defined(HINT)
		cout << "Query Times = " << query_times << endl;
		cout << "Query Avg Time = " << (double)(sumTime / (double)query_times) / 1000 << "ms" << endl;
		cout << endl << endl;
#endif
	}
    //至此每一个测试文件都遍历执行结束
	delete[] dataArray;
	return 0;
}
int edgeExistenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）
    string input_file_prefix = "";
	string input_file_suffix = "";
	string output_file_prefix = "";
	string output_file_suffix = "";
	string time_file_suffix = "_time.txt";
	switch (flag)//flag的作用仅是区别文件名，调用的方法是一样的
	{
	case 1:
		input_file_prefix = "_random_edge_existence_";
		input_file_suffix = "_sorted.txt";
		output_file_prefix = "_edge_existence_baseline_";
		output_file_suffix = "_res.txt";
		break;
	case 2:
		input_file_prefix = "_bool_";
		input_file_suffix = ".txt";
		output_file_prefix = "_bool_baseline_";
		output_file_suffix = "_res.txt";
		break;
	default:
		break;
	}

	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
    //2.一个一个测试文件按顺序遍历、读取数据项、执行边权重查询
	for (int i = 0; i < num.size(); i++) {
		//2.1读取当前测试文件内容到dataArray
		int datanum = readRandomFileToDataArray(input_dir + dataset_name + input_file_prefix + to_string(num[i]) + input_file_suffix, dataArray);

#if defined(DEBUG) || defined(HINT)
		cout << "****************** timeslot = " << num[i] << " ******************" << endl;
#endif
        //2.2创建当前测试文件对应的输出文件
		ofstream resultFile, timeFile;
		if (write) {
			char dir_path[FILENAME_MAX];
			strcpy(dir_path, output_dir.c_str());
			if (createDirectory(dir_path) != 0) {
				cout << "createDirectory error" << endl;
				return -1;
			}
			resultFile.open(output_dir + dataset_name + output_file_prefix + to_string(num[i]) + output_file_suffix);
			if (!resultFile.is_open()) {
				cout << "error in open file " << (output_dir + dataset_name + output_file_prefix + to_string(num[i]) + output_file_suffix) << endl;
				return -1;
			}
			timeFile.open(output_dir + dataset_name + output_file_prefix + to_string(num[i]) + time_file_suffix);
			if (!timeFile.is_open()) {
				cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + to_string(num[i]) + time_file_suffix) << endl;
				return -1;
			}
		}

        //2.3开始对当前dataArray中的每一个查询项执行边权重查询
		double sumTime = 0, sumTime_perquery = 0;
		int ones = 0;
		timeval tp1, tp2;
		for (int m = 0; m < query_times; m++) {
			sumTime_perquery = 0;
			for (int n = 0; n < datanum; n++) {
				gettimeofday( &tp1, NULL);
				//对文件中的每一个查询项调用edgeFrequenceBaseline方法
				int64_t res = edgeFrequenceBaseline(*dstgs, dataArray[n].source, dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
				gettimeofday( &tp2, NULL);
				double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
				sumTime_perquery += delta_t;
				if (write && (m == 0)) {
					if (res > 0)   
						ones++;
					if(n == (datanum - 1)) {
						resultFile << ((res > 0) ? 1 : 0);
						timeFile  << delta_t;
						break;
					}
					else {
						resultFile << ((res > 0) ? 1 : 0) << endl;
						timeFile  << delta_t << endl;
					}
				}
			}
			sumTime += (sumTime_perquery / (double)datanum);
		}
		if (write) {
			resultFile.flush();
			timeFile.flush();
			resultFile.close();
			timeFile.close();
		}
#if defined(DEBUG) || defined(HINT)
		cout << "ones: " << ones << endl;
		cout << "Query Times = " << query_times << endl;
		cout << "Query Avg Time = " << (double)(sumTime / (double)query_times) / 1000 << "ms" << endl;
		cout << endl << endl;
#endif
	}
    //至此每一个测试文件都遍历执行结束
	delete[] dataArray;
	return 0;
}
int nodeFrequenceBaselineTest_seq(string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line){
    //1.根据当前数据集，确定输入输出文件名的前后缀（为读取和写出文件做准备）
    string input_file_prefix = "";
	string input_file_suffix = "";
	string output_file_prefix = "";
	string output_file_suffix = "";
	string time_file_suffix = "_time.txt";
	switch (flag)
	{
	case 1:
		input_file_prefix = "_random_node_frequence_in_";
		input_file_suffix = "_sorted.txt";
		output_file_prefix = "_node_frequence_in_";
		output_file_suffix = "_res.txt";
		break;
	case 2:
		input_file_prefix = "_random_node_frequence_out_";
		input_file_suffix = "_sorted.txt";
		output_file_prefix = "_node_frequence_out_";
		output_file_suffix = "_res.txt";
		break;
	default:
		break;
	}
	QueryPairData* dataArray = new QueryPairData[query_data_pairs];
    //2.一个一个测试文件按顺序遍历、读取数据项、执行边权重查询
	for (int i = 0; i < num.size(); i++) {
		//2.1读取当前测试文件内容到dataArray
		int datanum = readRandomFileToDataArray(input_dir + dataset_name + input_file_prefix + to_string(num[i]) + input_file_suffix,dataArray);
#if defined(DEBUG) || defined(HINT)
		cout << "****************** timeslot = " << num[i] << " ******************" << endl;
#endif

        //2.2创建当前测试文件对应的输出文件
		ofstream resultFile, timeFile;
		if (write) {
			char dir_path[FILENAME_MAX];
			strcpy(dir_path, output_dir.c_str());
			if (createDirectory(dir_path) != 0) {
				cout << "createDirectory error" << endl;
				return -1;
			}

			if(line == 0) {
				resultFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + output_file_suffix);
				if (!resultFile.is_open()) {
					cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + output_file_suffix) << endl;
					return -1;
				}
				timeFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + time_file_suffix);
				if (!timeFile.is_open()) {
					cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + time_file_suffix) << endl;
					return -1;
				}
			}
			else {	// append
				resultFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + output_file_suffix, ios::app);
				cout << "append " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + output_file_suffix) << endl;
				if (!resultFile.is_open()) {
					cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + output_file_suffix) << endl;
					return -1;
				}
				timeFile.open(output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + time_file_suffix, ios::app);
				cout << "append " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + time_file_suffix) << endl;
				if (!timeFile.is_open()) {
					cout << "Error in open file, Path = " << (output_dir + dataset_name + output_file_prefix + "baseline_" + to_string(num[i]) + time_file_suffix) << endl;
					return -1;
				}
			}
		}

        //2.3开始对当前dataArray中的每一个查询项执行边权重查询
		double sumTime = 0, sumTime_perquery = 0;
		timeval tp1, tp2;
		for (int m = 0; m < query_times; m++) {
			sumTime_perquery = 0;
			for (int n = 0; n < datanum; n++) {
				if((m == 0) && (n < line)) 
					continue;
				gettimeofday( &tp1, NULL);
				//对每一个查询项调用nodeFrequenceBaseline方法
				int64_t res = nodeFrequenceBaseline(*dstgs, dataArray[n].source, (int)dataArray[n].destination, dataArray[n].start_time, dataArray[n].end_time);
				gettimeofday( &tp2, NULL);
				double delta_t = (tp2.tv_sec - tp1.tv_sec) * 1000000 +  (tp2.tv_usec - tp1.tv_usec);
				sumTime_perquery += delta_t;
				if (write && m == 0) {
					if(n == (datanum - 1)) {
						resultFile << res;
						timeFile  << delta_t;
						break;
					}
					else {
						resultFile << res << endl;
						timeFile  << delta_t << endl;
					}
				}
			}
			sumTime += (sumTime_perquery / (double)datanum);
		}
		if (write) {
			resultFile.flush();
			timeFile.flush();
			resultFile.close();
			timeFile.close();
		}
#if defined(DEBUG) || defined(HINT)
		cout << "Query Times = " << query_times << endl;
		cout << "Query Avg Time = " << (double)(sumTime / (double)query_times) / 1000 << "ms" << endl;
		cout << endl << endl;
#endif
	}
    //至此每一个测试文件都遍历执行结束
	delete[] dataArray;
	return 0;
}
//the functions that called by main
void edgeFrequenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write) {
	if(para_query)
		edgeFrequenceBaselineTest_para(input_dir, output_dir, dataset_name, num, query_times, write);
	else 
		edgeFrequenceBaselineTest_seq(input_dir, output_dir, dataset_name, num, query_times, write);
}
void edgeExistenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag) {
	if(para_query)
		edgeExistenceBaselineTest_para(input_dir, output_dir, dataset_name, num, query_times, write, flag);
	else 
		edgeExistenceBaselineTest_seq(input_dir, output_dir, dataset_name, num, query_times, write, flag);
}
void nodeFrequenceBaselineTest(bool para_query, string input_dir, string output_dir, string dataset_name, vector<int> num, int query_times, bool write, int flag, int line) {
	if(para_query)
		nodeFrequenceBaselineTest_para(input_dir, output_dir, dataset_name, num, query_times, write, flag, line);
	else 
		nodeFrequenceBaselineTest_seq(input_dir, output_dir, dataset_name, num, query_times, write, flag, line);
}
#endif