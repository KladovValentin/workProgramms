#include <iostream>
#include <vector>
#include <string>
#include<Windows.h>
using namespace std;

struct Constants {
	static const std::string s1;
	static const std::string s2;
};

struct knot {
public:
	knot(int k, int d, knot* l, knot* r) {
		key = k; data = d; knot* left_key = l; knot* right_key = r;
	}

	void SetRightKey(knot* r) {
		right_key = r;
	}
	knot GetRightKey() {
		return *right_key;
	}
	int pkey() {
		return key;
	}

private:
	int key = 1;
	int data = 1;
	knot* left_key = NULL;
	knot* right_key = NULL;
};

int findeq() {
	int n = 500;
	std::vector<int> nums;
	for (int i = 0; i < n + 1; i++) {
		nums.push_back(rand() % n + 1);
		std::cout << nums[i] << " ,";
	}
	std::cout << std::endl;
	int tortoise = nums[0];
	int hare = nums[0];
	while (1) {
		tortoise = nums[tortoise];
		hare = nums[nums[hare]];
		if (tortoise == hare)
			break;
	}
	int ptr1 = nums[0];
	int ptr2 = tortoise;
	while (ptr1 != ptr2) {
		ptr1 = nums[ptr1];
		ptr2 = nums[ptr2];
	}

	std::cout << nums[tortoise] << std::endl;
	return 0;
}

void testmap() {
	std::vector<knot> tree;
	tree.push_back({0,0,NULL,NULL});
	for (int i = 1; i < 10; i++) {
		tree.push_back({i,i*10,&tree[i-1],NULL});
	}
	for (int i = 0; i < 9; i++) {
		tree[i].SetRightKey(&tree[i + 1]) ;
	}
	for (int i = 0; i < 9; i++) {
		std::cout << tree[i].GetRightKey().pkey() << std::endl;
	}
}


int main() {
	cout << "|";
	for (int i = 0; i < 10000 / 500; i++)
		cout << " ";
	cout << "|" << endl;
	cout << " ";

	int Count = 0, count1 = 0;
	for (int e = 0; e < 10000; e++) {
		Count += 1;
		if (Count == 500) {
			count1 += 1;
			cout << "-";
			Count = 0;
			Sleep(1000000);
		}
	}
	cout << endl;
}