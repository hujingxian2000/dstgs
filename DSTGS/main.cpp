#include "QueryFunction.h"
#include <iomanip>
time_type getDatasetStartTime(string datasetPath);
int main(int argc, char* argv[]){
    cout<<fixed;
#if defined(DEBUG) || defined(HINT)  
	cout << setprecision(7);
	timeval main_start, main_end;
	gettimeofday( &main_start, NULL);

	// for (int i = 0; i < argc; i++) {
	// 	cout << argv[i] << " ";
	// }
	// cout << endl << endl;
#endif
    // 参数定义
	time_type startTime;						//数据集起始时间
	uint32_t width, depth;						//矩阵宽度 矩阵深度
	uint32_t granularityLength = 86400, gl = 1, slot = 2, fingerprintLen = 7;	//时间粒度长度，时间粒度，槽数，指纹长度
	string back_addr = "";

	int dataset = 1;							// 数据集编号
	int query_times = 1;						// 查询次数
	string filename, input_dir, output_dir;		// 数据集文件路径, 测试文件夹路径 , 输出文件夹路径
	string dataset_name;						// 数据集名称
	vector<int> num;							// num向量，表示L的集合
	int efflag = 0, eeflag = 0, nfflag = 0; 	// 边权重查询, 边存在查询, 点权重查询
	bool writeflag = false;						// 是否写入文件
	int node_query_flag = 0;					// 1-node_in_query, 2-node_out_query
	int edge_existence_flag = 1;				// 1-edge_existence_query, 2-bool_query
	int line = 0;								// 0-覆盖写入，1-追加写入
	bool para_query = true;  					// 0-sequential query, 1-parallel query

	uint32_t row_addrs = 4, column_addrs = 4;	//行地址数，列地址数
	bool edge_file_test = false;				
	bool node_file_test = false;

    // 命令行参数解析
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-dataset") == 0) {
			dataset = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-gl") == 0) {
			granularityLength = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-slot") == 0) {
			slot = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-fplength") == 0) {
			fingerprintLen = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-edgeweight") == 0) {
			efflag = 1;
		}
		if (strcmp(argv[i], "-edgeexistence") == 0) {
			eeflag = 1;
		}
		if (strcmp(argv[i], "-nodeinweight") == 0) {
			nfflag = 1;
			node_query_flag = 1;
		}
		if (strcmp(argv[i], "-nodeoutweight") == 0) {
			nfflag = 1;
			node_query_flag = 2;
		}
		if (strcmp(argv[i], "-bool") == 0) {
			edge_existence_flag = 2;
		}
		if (strcmp(argv[i], "-write") == 0) {
			writeflag = true;
		}
		if (strcmp(argv[i], "-para_query") == 0) {
			para_query = true;
		}
		if (strcmp(argv[i], "-seq_query") == 0) {
			para_query = false;
		}
		if (strcmp(argv[i], "-row_addrs") == 0) {
			row_addrs = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-col_addrs") == 0) {
			column_addrs = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-qtimes") == 0) {
			query_times = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-line") == 0) {		
			line = atoi(argv[++i]);
		}
	}
    switch (dataset) {
		case 1:
			filename = "..//Dataset//lkml";
			input_dir = "..//TestFiles//lkml//input//";
			output_dir = "..//TestFiles//lkml//output//";
			dataset_name = "lkml";
			num = { 8, 16, 32, 64, 128, 256, 512, 1024, 1536, 2048, 2560 };
			width = 2680;
			depth = 2680;
			break;
		case 2:
			filename = "..//Dataset//wiki-talk";
			input_dir = "..//TestFiles//wiki-talk//input//";
			output_dir = "..//TestFiles//wiki-talk//output//";
			dataset_name = "wiki-talk";
			num = { 32, 64, 128, 256, 512, 1024, 2048, 3072, 4096, 5120 };
			width = 13500;
			depth = 13500;
			break;
		case 3:
			filename = "..//Dataset//stackoverflow";
			input_dir = "..//TestFiles//stackoverflow//input//";
			output_dir = "..//TestFiles//stackoverflow//output//";
			dataset_name = "stackoverflow";
			num = { 8, 16, 32, 64, 128, 256, 512, 1024, 1536, 2048, 2560 };
			width = 20396;
			depth = 20396;
			break;
		case 4:
			filename = "..//Dataset//caida";
			input_dir = "..//TestFiles//caida//input//";
			output_dir = "..//TestFiles//caida//output//";
			dataset_name = "caida";
			num = { 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192, 9216, 10240, 11264, 12288 };
			width = 55000;
			depth = 55000;
			break;
		default:
			break;
	}
    for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-vector") == 0) {
			num.clear();
			while(true) {
				if((i + 1) >= argc) {
					break;
				}
				if(argv[i + 1][0] < '0' || argv[i + 1][0] > '9') {
					break;
				}
				else {
					num.push_back(atoi(argv[++i]));
				}
			}
		}
		if (strcmp(argv[i], "-input_dir") == 0) {
			input_dir = argv[++i];
			input_dir += "//";
		}
		if (strcmp(argv[i], "-output_dir") == 0) {
			output_dir = argv[++i];
			output_dir += "//";
		}
		if (strcmp(argv[i], "-filename") == 0) {
			filename = argv[++i];
		}
		if (strcmp(argv[i], "-w") == 0) {
			width = atoi(argv[++i]);
		}
		if (strcmp(argv[i], "-d") == 0) {
			depth = atoi(argv[++i]);
		}
	}


    // 获取当前数据集的起始时间
	startTime = getDatasetStartTime(filename);
	// 定义DSTGS_VAR类型的变量 dstgs_var
	DSTGS_VAR dstgs_var = { startTime, granularityLength, gl, width, depth, fingerprintLen, row_addrs, column_addrs };

#if defined(DEBUG) || defined(HINT)
	cout << "*******************************************************" << endl;
	cout << "Print Infomation" << endl;
	cout << "DSTGS_VAR: startTime = " << dstgs_var.startTime 
		 << ", granularityLength = " << dstgs_var.granularityLength 
		 << ", gl = " << dstgs_var.gl 
		 << ", width = " << dstgs_var.width 
		 << ", depth = " << dstgs_var.depth 
		 << ", fingerprintLen = " << dstgs_var.fingerprintLen << endl;
#endif

	//生成输出文件夹
	back_addr = "-DSTGS-" + to_string(row_addrs) + "x" + to_string(column_addrs) + "-res";
	string test_situation_dir = dataset_name + "_gl_" + to_string(dstgs_var.granularityLength) + "_" + to_string(dstgs_var.width) + "_" + 
		to_string(dstgs_var.depth) + "_" + to_string(SLOTNUM) + "_" + to_string(dstgs_var.fingerprintLen) + back_addr + "//";	
	output_dir += test_situation_dir;

	char dir_path[FILENAME_MAX];
	strcpy(dir_path, output_dir.c_str());
	if (createDirectory(dir_path) != 0) {
		cout << "Create Directory error" << endl;
		return -1;
	}

#if defined(DEBUG) || defined(HINT)
	cout << "dataset: " << filename << endl;
	cout << "input_dir: " << input_dir << endl;
	cout << "output_dir: " << output_dir << endl;
	cout << "write flag = " << writeflag << endl;
	cout << "vector num: ";
	for (int i = 0; i < num.size(); i++) {
		cout << num[i] << " ";
	}
	cout << endl;
	cout << "*******************************************************" << endl << endl;
#endif

	//正式开始流边处理：
	//边插入，调用baselineInsert方法
	cout << "****************** DSTGS insert start *****************" << endl;
	baselineInsert(dstgs_var, filename);	//传入DSTGS_VAR类型的变量和数据集路径
	cout << "****************** DSTGS insert end *******************" << endl << endl;
	
	//边权重查询，调用edgeFrequenceBaselineTest方法
	if (efflag == 1) {
		cout << "**************** DSTGS frequence start ****************" << endl;
		edgeFrequenceBaselineTest(para_query, input_dir, output_dir, dataset_name, num, query_times, writeflag);
		cout << "***************** DSTGS frequence end *****************" << endl << endl;
	}
	
	//边存在查询，调用edgeExistenceBaselineTest方法
	if (eeflag == 1) {
		cout << "**************** DSTGS existence start ****************" << endl;
		//最后一个参数代表查询类型（flag==1，代表存在性查询，flag==2代表可达性查询）
		edgeExistenceBaselineTest(para_query, input_dir, output_dir, dataset_name, num, query_times, writeflag, edge_existence_flag);
		cout << "***************** DSTGS existence end *****************" << endl << endl;
	}
	
	//点权重查询，调用nodeFrequenceBaselineTest方法
	if (nfflag == 1) {
		cout << "************* DSTGS node frequence start **************" << endl;
		//flag表示查询类型，line表示写入方式
		nodeFrequenceBaselineTest(para_query, input_dir, output_dir, dataset_name, num, query_times, writeflag, node_query_flag, line);
	}

	//输出整个程序执行时间
	gettimeofday( &main_end, NULL);
	double main_time = (main_end.tv_sec - main_start.tv_sec) + (main_end.tv_usec - main_start.tv_usec) / 1000000.0;
	cout << endl << "This program lasts for " << main_time / 60.0 << " min" << endl;
	return 0;

}
time_type getDatasetStartTime(string datasetPath) {
	ifstream ifs;
	ifs.open(datasetPath);
	if(!ifs.is_open()) {
		cout << "Open dataset error! Path = " << datasetPath << endl;
		return -1;
	}
	int64_t s, d;
	weight_type w;
	time_type startTime;
	ifs >> s >> d >> w >> startTime;
	ifs.close();
	if(startTime > 0) 
		return startTime - 1;
	else
		return -1;
}