#ifndef _DSTGS_H
#define _DSTGS_H
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <cmath>
#include <stdlib.h>
#include <malloc.h>
#include <bitset>
#include <memory.h>
#include <algorithm>
#include <sys/time.h>
#include "HashFunction.h"
#include "params.h"
using namespace std;


// 结构体FlatMap 
struct FlatMap {
    weight_type weights[FlatmapLen]; // 存储权重值
    uint32_t starttime;     // FlatMap 的起始时间
    uint32_t endtime;       // FlatMap 的结束时间

    FlatMap(uint32_t start, uint32_t end) : starttime(start), endtime(end) {
        for (int i = 0; i < FlatmapLen; ++i) {
            weights[i] = 0;
        }
    }

    // 插入
    void insert(int index, weight_type value) {
        if (index >= 0 && index < FlatmapLen) {
            weights[index] += value;
        }
    }

    // 查询
    weight_type get(int index){
        if (index >= 0 && index < FlatmapLen) {
            return weights[index];
        }
    }
};
//结构体basket，相当于矩阵当中的一个bucket
struct basket {
	uint16_t src[SLOTNUM];
	uint16_t dst[SLOTNUM]; 
    vector<FlatMap> flatMaps[SLOTNUM]; // 每个槽对应多个 FlatMap
    basket() {
        for (int i = 0; i < SLOTNUM; ++i) {
            flatMaps[i].emplace_back(1, CoverageLen); // 初始化第一个 FlatMap，覆盖时间范围为 [1, CoverageLen]
        }
    }
};
//结构体node，用于缓冲区，代表一个节点
struct node {
	uint32_t key;
	vector<FlatMap> flatMaps;
    node(){
        flatMaps.emplace_back(1, CoverageLen); // 初始化第一个 FlatMap，覆盖时间范围为 [1, CoverageLen]
    }
};
//结构体interval，子区间
struct interval {
    int64_t start;
    int64_t end;
};

// findx 和 findv 的核心作用是使得 findx 和 findv 类的实例能够作为函数对象（functor）使用，通常用于 find_if 等 STL 算法中查找符合条件的元素。
// findx 类用于查找单个 node，通过比较 node 的 key 是否与 value 匹配。
// findv 类用于查找包含多个 node 的 vector，通过比较 vector 中第一个 node 的 key 是否与 value 匹配。
class findx {
public:
	findx(const uint32_t va) { value = va; }
	bool operator()(const vector<node>::value_type &nod) {
		if (nod.key == value)
			return true;
		else
			return false;
	}
private:
	uint32_t value;
};
class findv {
public:
	findv(uint32_t va) { value = va; }
	bool operator()(vector<node> &vc) {
		if (vc[0].key == value)
			return true;
		else
			return false;
	}
private:
	uint32_t value;
};
class DSTGS
{
private:
	size_t matrixMemorySize;  // 用于记录矩阵占用的内存大小
    size_t bufferMemorySize;  // 用于记录邻接列表缓冲区占用的内存大小
    const uint32_t granularity;			// the granularity of the level
	uint32_t width;						// the width of the matrix
	uint32_t depth;						// the depth of the matrix
	const uint32_t fingerprintLength;	// the fingerprint length
	const uint32_t row_addrs;			// the row addrs
	const uint32_t column_addrs;		// the column_addrs
	basket* value;						// the matrix，value 是一个指向 basket 结构体数组的指针，表示图的矩阵。每个 basket 结构体代表一个图的矩阵位置（bucket），存储了多个源节点、目标节点以及对应的FlatMap。
	vector<vector<node>> successorAdjacencyList;//邻接列表缓冲区，是一个双重向量。
	int getMinIndex(uint32_t* a, int length);
public:
    DSTGS(uint32_t granularity, uint32_t width, uint32_t depth, uint32_t fingerprintLength, uint32_t row_addrs, uint32_t column_addrs);
    ~DSTGS();
	uint32_t getGranularity() const;
	void bucketCounting();
    void insert(string src, string dst, weight_type weight, uint32_t time);
    void flash();
	size_t getTotalMemorySize() const;  // 用于获取总的内存使用量（对外接口）
    void updateMemorySizes();  // 用于更新 matrixMemorySize 和 bufferMemorySize 的方法
	vector<interval> decomposeTargetRange(int64_t Tb, int64_t Te);
    weight_type edgeQuery(string src, string dst, uint32_t tb, uint32_t te);
    weight_type nodeQuery(string vertex, int type, uint32_t tb, uint32_t te);
    //bool reachabilityQuery(string s, string d);
};

// 构造函数
DSTGS::DSTGS(uint32_t granularity, uint32_t width, uint32_t depth, uint32_t fingerprintLength, uint32_t row_addrs, uint32_t column_addrs):
granularity(granularity), width(width), depth(depth), fingerprintLength(fingerprintLength), row_addrs(row_addrs), column_addrs(column_addrs) {
    //打印构造函数调用时的参数值
	cout << "DSTGS::DSTGS(granularity: " << granularity 
		<< ", width: " << width <<", depth: " << depth << ", fplen: " << fingerprintLength
		<< ", row_addrs: " << row_addrs << ", column_addrs: " << column_addrs << ")" << endl;
	//计算需要分配的内存大小 msize，用于存储如 basket 类型的数组
	uint32_t msize = width * depth;
	//按照 64 字节对齐，分配足够存储 msize 个 basket 结构体的内存空间，并将其赋值给 this->value。
	this->value = (basket *) memalign(64, sizeof(basket) * msize);

	// 检查内存分配是否成功，如果失败，输出错误信息并终止程序，避免后续出现空指针解引用等问题
    if (this->value == NULL) {
        cout << "Memory allocation failed in DSTGS constructor!" << endl;
        exit(EXIT_FAILURE);
    }
	cout << "内存分配成功" <<endl;

    // 初始化 value 所指向的每个 basket 结构体中的成员，采用更简单直接的方式
    for (uint32_t i = 0; i < msize; ++i) {
        // 对 basket 中的 src 和 dst 数组进行初始化，简单设为0
        std::fill_n(value[i].src, SLOTNUM, 0);
        std::fill_n(value[i].dst, SLOTNUM, 0);

        // 对于每个槽位对应的 flatMaps 向量，只进行简单的默认初始化，不做复杂操作
        for (int j = 0; j < SLOTNUM; ++j) {
            value[i].flatMaps[j] = std::vector<FlatMap>();
            // 如果需要确保至少有一个默认的 FlatMap，可以添加以下语句（这里先注释掉，可根据实际需求决定是否添加）
            // value[i].flatMaps[j].emplace_back(1, CoverageLen);
        }
    }
	//初始化
	//memset(this->value, 0, sizeof(basket) * msize);
	cout << "初始化成功" <<endl;
}
//析构函数
DSTGS::~DSTGS()
{
    cout << "DSTGS::~DSTGS()" << endl;
	free(this->value);
	vector<vector<node>>().swap(successorAdjacencyList);
}
int DSTGS::getMinIndex(uint32_t *a, int length)
{
    uint32_t min = a[0];
	int index = 0;
	for(int i = 1; i < length; i++) {
		if(a[i] < min) {
			min = a[i];
			index = i;
		}
	}
	return index;
}
uint32_t DSTGS::getGranularity() const
{
    if(this == NULL) {
		cout << "NULL pointer!!" << endl;
		getchar();
		exit(-1);
	}
	return this->granularity;
}
void DSTGS::bucketCounting()
{
	cout << "---------------------------------------" << endl;
	cout << "DSTGS bucketCounting(): print bucket..." << endl;
	int64_t room_count = 0;
	int64_t bucket_count = 0;
	for (int64_t i = 0; i < width * depth; i++) {
		if ((value[i].src[0] != 0) && (value[i].flatMaps[0][0].weights[0] != 0)) {
				bucket_count++;
		}
		for (int64_t j = 0; j < SLOTNUM; j++) {
			if ((value[i].src[j] != 0) && (value[i].flatMaps[j][0].weights[0] != 0)) {
				room_count++;
			}
		}
	}
	cout << "DSTGS room_count = " << room_count << ", total room = " << (width * depth * SLOTNUM) << ", space usage is " << 
			(double)room_count / (double)(width * depth * SLOTNUM) * 100 << "%" << endl;
	cout << "DSTGS bucket_count = " << bucket_count << ", total bucket = " << (width * depth) << ", space usage is " << 
			(double)bucket_count / (double)(width * depth) * 100 << "%" << endl;
	
	//print buffer size
	cout << "---------------------------------------" << endl;
	cout << "print successorAdjacencyList..." << endl;
	cout << "successorAdjacencyList.size() = " << successorAdjacencyList.size() << endl;
	cout << "successorAdjacencyList.capacity() = " << successorAdjacencyList.capacity() << endl;

	int64_t total_sucBuffer = 0;
	int64_t total_cap = 0, total_size = 0;
	for(int64_t i = 0; i < this->successorAdjacencyList.size(); i++) {
		total_cap += this->successorAdjacencyList[i].capacity();
		total_size += this->successorAdjacencyList[i].size();
	}
	cout << "total_size = " << total_size << endl;
	cout << "total_cap = " << total_cap << endl;
	for(int64_t i = 0; i < this->successorAdjacencyList.size(); i++) {
		 total_sucBuffer += this->successorAdjacencyList[i].size();
	}
	cout << "total_sucBuffer = " << total_sucBuffer << endl;
	cout << "---------------------------------------" << endl;
	return;
}
void DSTGS::updateMemorySizes()
{
    // 计算矩阵（value 指向的内存区域）占用的内存大小
    matrixMemorySize = width * depth * sizeof(basket);
    // 计算邻接列表缓冲区（successorAdjacencyList）占用的内存大小
    bufferMemorySize = 0;
    for (const auto& subList : successorAdjacencyList)
    {
        bufferMemorySize += subList.size() * sizeof(node);
        for (const auto& nodeInstance : subList)
        {
            bufferMemorySize += nodeInstance.flatMaps.size() * sizeof(FlatMap);
            // 还可以考虑 FlatMap 内部成员数组 weights 的内存占用情况，如果需要更精确统计的话
             bufferMemorySize += nodeInstance.flatMaps.size() * FlatmapLen * sizeof(weight_type);
        }
    }
}
size_t DSTGS::getTotalMemorySize() const
{
    return matrixMemorySize + bufferMemorySize + sizeof(granularity) + sizeof(width) + sizeof(depth) + sizeof(fingerprintLength)+ sizeof(row_addrs) + sizeof(column_addrs);
}

// 边插入方法
void DSTGS::insert(string src, string dst, weight_type weight, uint32_t time)
{
    // 使用哈希函数 hfunc[HASH] 计算源节点和目标节点的哈希值——相当于H(s)和H(v)：
    uint32_t hash_src = (*hfunc[HASH])((unsigned char *)(src.c_str()), src.length());
    uint32_t hash_dst = (*hfunc[HASH])((unsigned char *)(dst.c_str()), dst.length());
    
    // 下面计算源节点和目的节点的指纹——相当于f(s)和f(v)：
    uint32_t mask = (1 << fingerprintLength) - 1;
    uint16_t fp_src = hash_src & mask;
    if (fp_src == 0) fp_src += 1;
    uint16_t fp_dst = hash_dst & mask;
    if (fp_dst == 0) fp_dst += 1;
    
    // 下面计算源节点和目的节点的哈希地址——相当于h(s)和h(v)：
    uint32_t addr_src = (hash_src >> fingerprintLength) % depth;
    uint32_t addr_dst = (hash_dst >> fingerprintLength) % width;
    
    //k1 和 k2 是将矩阵地址和指纹组合后生成的唯一键值，用于快速标识源节点和目标节点。（k1和k2都是在缓冲区中会利用的）
	uint32_t k1 = (addr_src << fingerprintLength) + fp_src;
	uint32_t k2 = (addr_dst << fingerprintLength) + fp_dst;
    
    // 下面进行种子的生成：
	uint32_t head = 16384;								// pow(2, 14);
	uint32_t* seed1 = new uint32_t[row_addrs];			// seed1：行地址种子数组。长度为 row_addrs。
	uint32_t* seed2 = new uint32_t[column_addrs];		// seed2：列地址种子数组。长度为 column_addrs。
	// fp_src 和 fp_dst分别为源节点和目标节点的指纹，用作种子数组的初始值。
	seed1[0] = fp_src;
	seed2[0] = fp_dst;
	// 经典的线性同余生成
	for (int i = 1; i < row_addrs; i++)	
		seed1[i] =  (seed1[i - 1] * multiplier + increment) % modulus;
	for (int i = 1; i < column_addrs; i++)	
		seed2[i] =  (seed2[i - 1] * multiplier + increment) % modulus;


	// 开始进行边插入——插入矩阵
	for (int i = 0; i < row_addrs; i++) {
		// 生成的 seed1 和 seed2 将用于计算哈希地址序列，得到row_addr（个数为row_addrs）和column_addr（个数为column_addrs）
		uint32_t row_addr = (addr_src + seed1[i]) % depth;
		for (int j = 0; j < column_addrs; j++) {
			uint32_t column_addr = (addr_dst + seed2[j]) % width;
			// pos 就是实际定位到的位置（矩阵中的一维索引）
			uint32_t pos = row_addr * width + column_addr;

			//遍历矩阵pos位置的每个插槽
			for (int m = 0; m < SLOTNUM; m++) {
				//如果索引匹配、指纹匹配，插入当前边
				if (((value[pos].src[m] >> 14) == i) && ((value[pos].dst[m] >> 14) == j) && ((value[pos].src[m] & mask) == fp_src) && ((value[pos].dst[m] & mask) == fp_dst)) {
                    // 先定位到当前pos当前槽的flatMaps(一组)
                    auto &flatMaps = value[pos].flatMaps[m];
                    // 计算 flatMaps 的索引号 k
                    int k = (time - 1) / CoverageLen;
                    // 判断是否需要新建 FlatMap
                    while (k >= flatMaps.size()) {
                        uint32_t newStartTime = flatMaps.back().endtime + 1;
                        uint32_t newEndTime = newStartTime + CoverageLen - 1;
                        flatMaps.emplace_back(newStartTime, newEndTime);
                    }
                    // 计算元素下标 ii
                    int ii = time + (1-k) * CoverageLen - 2;
                    // 更新对应的 FlatMap
                    flatMaps[k].insert(ii, weight);                   
                    
                    delete[] seed1;
					delete[] seed2;
					return;
				}
				//如果找到空槽，初始化后，插入当前边
				if (value[pos].src[m] == 0 && value[pos].dst[m] == 0) {
					value[pos].src[m] = fp_src + head * i;
					value[pos].dst[m] = fp_dst + head * j;

                    // 先定位到当前pos当前槽的flatMaps(一组)
                    auto &flatMaps = value[pos].flatMaps[m];
                    // 计算 flatMaps 的索引号 k
                    int k = (time - 1) / CoverageLen;
                    // 判断是否需要新建 FlatMap
                    while (k >= flatMaps.size()) {
						// 如果 flatMaps 为空，需要先初始化一个默认的 FlatMap，设置合适的起始时间和结束时间
						if (flatMaps.empty()) {
							flatMaps.emplace_back(1, CoverageLen);
						} else {
							uint32_t newStartTime = flatMaps.back().endtime + 1;
							uint32_t newEndTime = newStartTime + CoverageLen - 1;
							flatMaps.emplace_back(newStartTime, newEndTime);
						}
						// cout << "新建FlatMap" <<endl;
                        // uint32_t newStartTime = flatMaps.back().endtime + 1;
                        // uint32_t newEndTime = newStartTime + CoverageLen - 1;
                        // flatMaps.emplace_back(newStartTime, newEndTime);
                    }
                    // 计算元素下标 ii
                    int ii = time + (1-k) * CoverageLen - 2;
                    // 更新对应的 FlatMap
                    flatMaps[k].insert(ii, weight); 

					delete[] seed1;
					delete[] seed2;
					return;
				}
			}
		}
	}
	delete[] seed1;
	delete[] seed2;

    // 如果矩阵插入失败，将当前 k1 和 k2 的关系保存到缓冲区 successorAdjacencyList。
	// successorAdjacencyList 是一个 vector<vector<node>>，表示邻接表的形式。每个 vector<node> 存储一个节点及其相关的后续节点。node 的结构包含：key：保存节点标识。FlatMaps：存储权重。
	// 首先，使用 find_if 方法在 successorAdjacencyList 中寻找 k1：
	vector<vector<node> >::iterator it = find_if(successorAdjacencyList.begin(), successorAdjacencyList.end(), findv(k1));
	// 如果找到了k1
	if (it != successorAdjacencyList.end()) {		
		// 继续在邻接列表中查找是否已有 k2
		vector<node>::iterator iter = find_if(it->begin(), it->end(), findx(k2));
		// 如果找到 k2：更新权重 weight。
		if (iter != it->end()) {	
			//iter->weight += weight;
            // 先定位到当前node的flatMaps(一组)
            auto &flatMaps=iter->flatMaps;
            // 计算 flatMaps 的索引号 k
            int k = (time - 1) / CoverageLen;
            // 判断是否需要新建 FlatMap
            while (k >= flatMaps.size()) {
                uint32_t newStartTime = flatMaps.back().endtime + 1;
                uint32_t newEndTime = newStartTime + CoverageLen - 1;
                flatMaps.emplace_back(newStartTime, newEndTime);
            }
            // 计算元素下标 ii
            int ii = time + (1-k) * CoverageLen - 2;
            // 更新对应的 FlatMap
            flatMaps[k].insert(ii, weight); 
		}
		// 如果未找到 k2：创建一个新节点 tmpnode，保存 k2 和flatMaps，添加到 k1 的邻接列表中。
		else {
			node tmpnode;
			tmpnode.key = k2;
			//tmpnode.weight = weight;
            // 先定位到当前node的flatMaps(一组)
            auto &flatMaps=tmpnode.flatMaps;
            // 计算 flatMaps 的索引号 k
            int k = (time - 1) / CoverageLen;
            // 判断是否需要新建 FlatMap
            while (k >= flatMaps.size()) {
                uint32_t newStartTime = flatMaps.back().endtime + 1;
                uint32_t newEndTime = newStartTime + CoverageLen - 1;
                flatMaps.emplace_back(newStartTime, newEndTime);
            }
            // 计算元素下标 ii
            int ii = time + (1-k) * CoverageLen - 2;
            // 更新对应的 FlatMap
            flatMaps[k].insert(ii, weight); 

			it->push_back(tmpnode);
		}
	}
    //如果没有找到k1
	else {
		// 新建邻接列表：创建一个 vector<node>，将 k1 作为首节点。
		node tmpnode;
		tmpnode.key = k1;
		//tmpnode.weight = 0;
		vector<node> vc;
		vc.push_back(tmpnode);

		//检查是否存在自环：
        //如果 k1 != k2（不是自环），添加 k2 节点。
		if (k1 != k2) {
			node newnode;
			newnode.key = k2;
			//newnode.weight = weight;
            // 先定位到当前node的flatMaps(一组)
            auto &flatMaps=newnode.flatMaps;
            // 计算 flatMaps 的索引号 k
            int k = (time - 1) / CoverageLen;
            // 判断是否需要新建 FlatMap
            while (k >= flatMaps.size()) {
                uint32_t newStartTime = flatMaps.back().endtime + 1;
                uint32_t newEndTime = newStartTime + CoverageLen - 1;
                flatMaps.emplace_back(newStartTime, newEndTime);
            }
            // 计算元素下标 ii
            int ii = time + (1-k) * CoverageLen - 2;
            // 更新对应的 FlatMap
            flatMaps[k].insert(ii, weight);  
			
            vc.push_back(newnode);
		}
        //如果 k1 == k2（出现自环），直接增加权重。
		else {
			//vc[0].weight += weight;
            // 先定位到当前node的flatMaps(一组)
            auto &flatMaps=vc[0].flatMaps;
            // 计算 flatMaps 的索引号 k
            int k = (time - 1) / CoverageLen;
            // 判断是否需要新建 FlatMap
             while (k >= flatMaps.size()) {
                uint32_t newStartTime = flatMaps.back().endtime + 1;
                uint32_t newEndTime = newStartTime + CoverageLen - 1;
                flatMaps.emplace_back(newStartTime, newEndTime);
            }
            // 计算元素下标 ii
            int ii = time + (1-k) * CoverageLen - 2;
            // 更新对应的 FlatMap
            flatMaps[k].insert(ii, weight);  
		}
		//将新的邻接列表 vc 添加到 successorAdjacencyList。
		successorAdjacencyList.push_back(vc);
	}
	return;

}
// Flash方法
void DSTGS::flash()
{
    //遍历矩阵的每一个pos
    for (uint32_t pos = 0; pos < width * depth; ++pos)
    {
        //遍历当前pos中的每一个槽
        for (int slot = 0; slot < SLOTNUM; ++slot)
        {
            //定位当前FlatMap（一组）
            vector<FlatMap>& flatMaps = value[pos].flatMaps[slot];
            //遍历每一个单独的FlatMap
            for(int k=0;k<flatMaps.size();k++)
            {
                for (int level = CoverageLen / 2; level > 0; level /= 2)//外层循环定层级
                {
                    int x=0;
                    for (int i = level - 1; i < 2 *level - 1; i++)//内层循环定下标
                    {                       
                        flatMaps[k].weights[i]=flatMaps[k].weights[i+level+x]+flatMaps[k].weights[i+level+x+1];//更新（两者求和）
                        x++;
                    }
                }
                
            }
        }
    }
}
//目标范围分解方法
vector<interval> DSTGS::decomposeTargetRange(int64_t Tb, int64_t Te) {
    vector<interval> subinterval;

    // 第一步，按FlatMap结构拆分目标范围，优化重复计算部分
    int64_t currentTb = Tb;
    int64_t currentTe = Te;
    vector<interval> subrange;
    while (true) {
        int64_t queryRangeLength = currentTe - currentTb + 1;
        bool needSplit = (queryRangeLength > CoverageLen) ||
                         (queryRangeLength < CoverageLen && currentTe > ((currentTb + CoverageLen - 1) / CoverageLen) * CoverageLen);

        if (needSplit) {
            int64_t r = ((currentTb + CoverageLen - 1) / CoverageLen) * CoverageLen;

            interval sub1 = {currentTb, r};
            subrange.push_back(sub1);

            currentTb = r + 1;
        } else {
            interval sub = {currentTb, currentTe};
            subrange.push_back(sub);
            break;
        }
    }

    // 第二步，针对每个初步拆分后的子范围进一步细分，优化条件判断和循环逻辑，重点调整步骤b及后续处理逻辑
    for (const auto& sub : subrange) {
        int64_t tb = sub.start;
        int64_t te = sub.end;
        int64_t l = te - tb + 1;

        // 利用位运算优化步骤a判断
        if (l == 1) {
            subinterval.push_back({tb, te});
            continue;
        }

        // 精准处理步骤b逻辑，针对偶数起始位置进行正确拆分
        if (tb % 2 == 0) {
            interval newSub1 = {tb, tb};
            subinterval.push_back(newSub1);
            tb++;
            l--;
            if (l == 1) {
                interval newSub2 = {tb, te};
                subinterval.push_back(newSub2);
                continue;
            }
        }

        // 步骤c逻辑基本保持，结合位运算简化判断，判断是否满足特殊情况直接拆分整个范围
        int64_t n = static_cast<int64_t>(log2(l));
        if (l == (1 << n) && (tb % l) == 1) {
            subinterval.push_back({tb, te});
            continue;
        }

        // 步骤d逻辑优化，合理拆分剩余范围
        while (l > 1) {
            int64_t maxLength = 1 << (63 - __builtin_clzll(l));
            while (maxLength > 1 && (tb % maxLength)!= 1) {
                maxLength >>= 1;
            }
            interval newSub1 = {tb, tb + maxLength - 1};
            subinterval.push_back(newSub1);

            tb += maxLength;
            l = te - tb + 1;
        }
    }
    return subinterval;
}
// 边查询方法
weight_type DSTGS::edgeQuery(string src, string dst, uint32_t start, uint32_t end)
{
	weight_type edgefrequent=0;
    //生成哈希值
	uint32_t hash_src = (*hfunc[HASH])((unsigned char*)(src.c_str()), src.length());
	uint32_t hash_dst = (*hfunc[HASH])((unsigned char*)(dst.c_str()), dst.length());
	
    //生成指纹和哈希地址
	uint32_t mask = (1 << fingerprintLength) - 1;
	uint16_t fp_src = hash_src & mask;
	if (fp_src == 0) fp_src += 1;
	uint32_t addr_src = (hash_src >> fingerprintLength) % depth;
	uint16_t fp_dst = hash_dst & mask;
	if (fp_dst == 0) fp_dst += 1;
	uint32_t addr_dst = (hash_dst >> fingerprintLength) % width;
	
    //生成种子
	uint32_t* seed1 = new uint32_t[row_addrs];			// row address seeds
	uint32_t* seed2 = new uint32_t[column_addrs];		// column address seeds
	seed1[0] = fp_src;
	seed2[0] = fp_dst;
	for (int i = 1; i < row_addrs; i++)	
		seed1[i] =  (seed1[i - 1] * multiplier + increment) % modulus;
	for (int i = 1; i < column_addrs; i++)	
		seed2[i] =  (seed2[i - 1] * multiplier + increment) % modulus;

	//开始查询，首先在矩阵中进行查询
	for (int i = 0; i < row_addrs; i++) {
		//首先得到row_addr（个数为row_addrs）和column_addr（个数为column_addrs）
		uint32_t row_addr = (addr_src + seed1[i]) % depth;
		for (int j = 0; j < column_addrs; j++) {
			uint32_t column_addr = (addr_dst + seed2[j]) % width;
			//然后定位到具体位置 pos
			uint32_t pos = row_addr * width + column_addr;
			if(pos >= width * depth || pos < 0) {
				cout << "matrix pos: " << pos << " out of range!" << endl;
				continue;
			}

			//对该pos中的每一个槽进行遍历
			for (int m = 0; m < SLOTNUM; m++) {
				//如果索引和指纹都匹配，那么找到该边，返回对应的权重值
				if (((value[pos].src[m] >> 14) == i) && ((value[pos].dst[m] >> 14) == j) && ((value[pos].src[m] & mask) == fp_src) && ((value[pos].dst[m] & mask) == fp_dst)) {
					// 先定位到当前pos当前槽的flatMaps(一组)
                    auto &flatMaps = value[pos].flatMaps[m];

					// 接下来进行区间分解
					//vector<interval> subranges=decomposeTargetRange(start,end);
    				//vector<interval> subintervals=furtherDecompose(subranges);
					vector<interval> subintervals=decomposeTargetRange(start,end);
					//遍历每一个可查询子区间
					for (const auto& sub : subintervals) {
						//result += DSTGS.edgeQuery(s,d,sub.start,sub.end);
						int l = sub.end - sub.start + 1;
						// 计算 flatMaps 的索引号 k
						int k = (sub.start - 1) / CoverageLen;
						// 计算元素下标 ii                   
						int ii = CoverageLen / l + ceil((double)(sub.start - k * CoverageLen) / (double) l) - 2;
						edgefrequent += flatMaps[k].get(ii);
					}

                    // 返回对应的权重
                    delete[] seed1;
					delete[] seed2;
                    return edgefrequent;                  

				}
			}
		}
	}
	delete[] seed1;
	delete[] seed2;

    //如果矩阵中没找到，查询缓冲区
	//将行地址与指纹结合形成 k1，将列地址与指纹结合形成 k2。k1 和 k2 共同表示 src -> dst 的唯一关系。（k1和k2就是在缓冲区使用的东西）
	uint32_t k1 = (addr_src << fingerprintLength) + fp_src;
	uint32_t k2 = (addr_dst << fingerprintLength) + fp_dst;
	// map<uint32_t, uint32_t>::iterator it = successorIndex.find(k1);
	//在 successorAdjacencyList 中查找 k1。
	vector<vector<node> >::iterator it = find_if(successorAdjacencyList.begin(), successorAdjacencyList.end(), findv(k1));
	//如果找到，继续查找 k2。
	if (it != successorAdjacencyList.end()) {
		vector<node>::iterator iter = find_if(it->begin(), it->end(), findx(k2));
		//如果找到匹配的 k2，返回对应权重。
		if (iter != it->end()) {
            // 先定位到当前的flatMaps(一组)
            auto &flatMaps = iter->flatMaps;
            // // 计算 flatMaps 的索引号 k
            // int k = (tb - 1) / CoverageLen;
            // // 计算元素下标 ii                   
            // int ii = CoverageLen / l + ceil((double)(tb - k * CoverageLen) / (double) l) - 2;
            // // 返回对应的权重
            // return flatMaps[k].get(ii);

			// 接下来进行区间分解
			vector<interval> subintervals=decomposeTargetRange(start,end);
			//遍历每一个可查询子区间
			for (const auto& sub : subintervals) {
				//result += DSTGS.edgeQuery(s,d,sub.start,sub.end);
				int l = sub.end - sub.start + 1;
				// 计算 flatMaps 的索引号 k
				int k = (sub.start - 1) / CoverageLen;
				// 计算元素下标 ii                   
				int ii = CoverageLen / l + ceil((double)(sub.start - k * CoverageLen) / (double) l) - 2;
				edgefrequent += flatMaps[k].get(ii);
			} 
			return edgefrequent; 
		}
	}
	//如果都没找到，返回0
	return 0;

}
// 点查询方法
weight_type DSTGS::nodeQuery(string vertex, int type, uint32_t tb, uint32_t te)
{
    //返回值weight
	weight_type weight = 0;
	//计算节点哈希值
	uint32_t hash_vertex = (*hfunc[HASH])((unsigned char*)(vertex.c_str()), vertex.length());
	//计算指纹
	uint32_t mask = pow(2, fingerprintLength) - 1;
	uint16_t fp = hash_vertex & mask;
	if (fp == 0) fp += 1;
	//如果type==0，则代表出度查询（定位行，遍历列），如果type==1，则代表入度查询（定位列，遍历行）
	int addrs = (type == 0) ? row_addrs : column_addrs;
	//计算种子
	uint32_t* seeds = new uint32_t[addrs];			// address seeds
	seeds[0] = fp;
	for (int i = 1; i < addrs; i++)	
		seeds[i] =  (seeds[i - 1] * multiplier + increment) % modulus;
	

	uint32_t l = te - tb + 1;
	//type==0，出度查询
	if (type == 0) {
		//首先查询矩阵
		//计算哈希地址
		uint32_t addr = (hash_vertex >> fingerprintLength) % depth;
		for (int i = 0; i < row_addrs; i++)	{
			//对每一个i生成哈希地址序列（行地址）
			uint32_t row_addr = (addr + seeds[i]) % depth;
			for (int j = 0; j < width; j++)	{
				//遍历所有列，生成具体位置pos
				uint32_t pos = row_addr * width + j;
				for (int m = 0; m < SLOTNUM; ++m) {
					//对该pos遍历所有槽，如果源节点的索引和指纹匹配，代表是同一个源节点
					if (((value[pos].src[m] >> 14) == i) && ((value[pos].src[m] & mask) == fp)) {
                        // 先定位到当前pos当前槽的flatMaps(一组)
                        auto &flatMaps = value[pos].flatMaps[m];
                        // 计算 flatMaps 的索引号 k
                        int k = (tb - 1) / CoverageLen;
                        // 计算元素下标 ii                   
                        int ii = CoverageLen / l + ceil((double)(tb - k * CoverageLen) / (double) l) - 2;
                        // 累积出度权重
						weight += flatMaps[k].get(ii);
					}
				}
			}	
		}

		//然后查询缓冲区
		//计算k1（源节点的唯一标识）
		uint32_t k1 = (addr << fingerprintLength) + fp;
		// map<uint32_t, uint32_t>::iterator it = successorIndex.find(k1);
		//在缓冲区中寻找k1
		vector<vector<node> >::iterator it = find_if(successorAdjacencyList.begin(), successorAdjacencyList.end(), findv(k1));
		if (it != successorAdjacencyList.end())	{
			//如果找到了k1，遍历所有记录，累加其权重。
			vector<node>::iterator iter;
			for (iter = it->begin(); iter != it->end(); iter++) {
				//weight += iter->weight;
                // 先定位到当前的flatMaps(一组)
                auto &flatMaps = iter->flatMaps;
                // 计算 flatMaps 的索引号 k
                int k = (tb - 1) / CoverageLen;
                // 计算元素下标 ii                   
                int ii = CoverageLen / l + ceil((double)(tb - k * CoverageLen) / (double) l) - 2;
                // 累积出度权重
				weight += flatMaps[k].get(ii);
			}
		}
	}
	//type==1，入度查询
	else if (type == 1) {
		//首先查询矩阵
		//计算哈希地址
		uint32_t addr = (hash_vertex >> fingerprintLength) % width;
		for (int i = 0; i < column_addrs; i++) {
			//对每一个i生成哈希地址序列（列地址）
			uint32_t col_addr = (addr + seeds[i]) % width;
			for (int j = 0; j < depth; j++) {
				//遍历所有行，生成具体位置pos
				uint32_t pos = j * width + col_addr;
				for (int m = 0; m < SLOTNUM; ++m) {
					//遍历该位置的所有槽，如果目的节点的索引和指纹匹配，代表是同一个目的节点
					if (((value[pos].dst[m] >> 14) == i) && ((value[pos].dst[m] & mask) == fp)) {
                        // 先定位到当前pos当前槽的flatMaps(一组)
                        auto &flatMaps = value[pos].flatMaps[m];
                        // 计算 flatMaps 的索引号 k
                        int k = (tb - 1) / CoverageLen;
                        // 计算元素下标 ii                   
                        int ii = CoverageLen / l + ceil((double)(tb - k * CoverageLen) / (double) l) - 2;
                        // 累积入度权重
						weight += flatMaps[k].get(ii);
					}
				}
			}
		}
		//然后查询缓冲区
		//计算k1，目标点的唯一标识（实际上是k2，只不过这里写k1而已）
		uint32_t k1 = (addr << fingerprintLength) + fp;
		//外层循环，遍历整个缓冲区
		for(vector<vector<node> >::iterator it = successorAdjacencyList.begin(); it != successorAdjacencyList.end(); it++){
			//内层循环，遍历当前邻接列表中的每一个 node，每个 node 存储了一个边的信息（即目标点和该边的权重）。
			for(vector<node>::iterator iter = it->begin(); iter!= it->end(); iter++){
				if(iter->key == k1){
					//检查当前 node 的 key 是否与查询的 k1 相同。如果相同，说明找到了与目标点 k1 对应的边。累积权重
					// 先定位到当前的flatMaps(一组)
                    auto &flatMaps = iter->flatMaps;
                    // 计算 flatMaps 的索引号 k
                    int k = (tb - 1) / CoverageLen;
                    // 计算元素下标 ii                   
                    int ii = CoverageLen / l + ceil((double)(tb - k * CoverageLen) / (double) l) - 2;
                    // 累积出度权重
                    weight += flatMaps[k].get(ii);
				}
			}
		}
	}
	delete[] seeds;
	return weight;
}

#endif