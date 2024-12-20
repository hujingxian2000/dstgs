#ifndef DECOMPOSEFUNCTION_H
#define DECOMPOSEFUNCTION_H
#include <iostream>
#include <vector>
#include <cmath>
#include "params.h"
using namespace std;

// 结构体
struct interval {
    int64_t start;
    int64_t end;
};

// 计算大于Tb的第一个FlatMap分界点
int64_t getFirstFlatMapBoundary(int64_t Tb) {
    return ((Tb + CoverageLen - 1) / CoverageLen) * CoverageLen;
}

// 分解目标查询范围
vector<interval> decomposeTargetRange(int64_t Tb, int64_t Te) {
    vector<interval> subrange;

    // 不断拆分直到每个子范围只涉及一个FlatMap结构
    int64_t currentTb = Tb;
    int64_t currentTe = Te;
    while (true) {
        int64_t queryRangeLength = currentTe - currentTb + 1;

        // 判断当前范围是否涉及多个FlatMap结构
        if ((queryRangeLength > CoverageLen) ||
            (queryRangeLength < CoverageLen && currentTe > getFirstFlatMapBoundary(currentTb))) {

            // 拆分当前范围
            int64_t r = getFirstFlatMapBoundary(currentTb);

            // 处理[currentTb, r]
            interval sub1;
            sub1.start = currentTb;
            sub1.end = r;
            subrange.push_back(sub1);

            // 更新当前范围为[r + 1, currentTe]，准备下一轮拆分（如果需要）
            currentTb = r + 1;
        } else {
            // 当前范围不涉及多个FlatMap结构，直接将其作为子范围添加到结果中
            interval sub;
            sub.start = currentTb;
            sub.end = currentTe;
            subrange.push_back(sub);
            break;
        }
    }
	//返回子范围集合subrange
    return subrange;
}

// 进一步分解
vector<interval> furtherDecompose(const vector<interval>& subrange) {
	vector<interval> subinterval;
	//遍历子范围集合中的每一个子范围
    for (const auto& sub : subrange) {
        int64_t tb = sub.start;
        int64_t te = sub.end;
        int64_t l = te - tb + 1;

        // 步骤a
		if (l == 1) {
            interval newSub = {tb, te};
            subinterval.push_back(newSub);
            continue;
        }
		// 步骤b
        if (tb % 2 == 0) {
            interval newSub = {tb, tb};
            subinterval.push_back(newSub);
            tb = tb + 1;
            l = te - tb + 1;
			// 跳转至步骤a
            if (l == 1) {
                interval newSub2 = {tb, te};
                subinterval.push_back(newSub2);
                continue;
            }
        }
		// 步骤c
        int64_t n = static_cast<int64_t>(log2(l));
        if (l == 1 << n && tb % l == 1) {
            interval newSub = {tb, te};
            subinterval.push_back(newSub);
            continue;
        }

		// 步骤d
        while (true) {
            int64_t maxLength = 1;
            while (maxLength <= l && (tb % maxLength == 1)) {
                maxLength *= 2;
            }
            maxLength /= 2;

            interval newSub1 = {tb, tb + maxLength - 1};
            subinterval.push_back(newSub1);

            if (tb + maxLength > te) {
                break;
            } else {
                tb = tb + maxLength;
                l = te - tb + 1;
            }
        }
    }
	return subinterval;
}
#endif