#ifndef _SHANGKE_TRIE_H_
#define _SHANGKE_TRIE_H_

#include <vector>
#include <string>

using namespace std;

class Trie {
	private:
		const static int CHAR_SET_NUM = 26;
		int cnt;
		int tot_cnt;
		Trie *next[CHAR_SET_NUM];
	public:
		/* *
		 * 构造器
		 */
		Trie();
		Trie(const Trie& o);

		/* *
		 * 析构函数
		 */
		~Trie();

		/* *
		 * 运算符重载
		 */
		// 合并两颗Trie树，返回新节点
		Trie operator + (const Trie& o);

		// 简单的复制，不增加新的空间
		Trie &operator = (const Trie& o);

		// 合并在原树上
		Trie &operator += (const Trie& o);

		/* *
		 * 一般函数
		 */
		// 整体拷贝复制，创建新节点
		void copy(const Trie& o);

		// 获取该单词的数目
		int get_cnt() const;

		// 获取该前缀的数目
		int get_tot_cnt() const;

		// 插入一个单词，返回单词词尾
		char* insert(char* s);

		// 查询某个单词的数目
		int query(char* s) const;

		// 获取某个前缀的数目
		int query_pre(char* s) const;

		// 删除某个单词，返回删除的数目
		int remove(char* s);

		// 获取所有的单词及其数目
		vector< pair<int, string> > get_all() const;

		// 将当前单词及其子孙插入到vector中
		void get_dfs(int deep, vector< pair<int, string> >& ans) const;
};

#endif
