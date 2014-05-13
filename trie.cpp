#include "trie.h"

#include <vector>
#include <string>
#include <algorithm>

using namespace std;

Trie::Trie() {
	cnt = 0, tot_cnt = 0;
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		next[i] = NULL;
	}
}

Trie::Trie(const Trie & o) {
	cnt = o.cnt, tot_cnt = o.tot_cnt;
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		next[i] = NULL;
	}
}

Trie::~Trie() {
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		if (next[i]) {
			delete next[i];
		}
	}
}

Trie Trie::operator + (const Trie& o) {
	Trie ans;
	ans.copy(*this), ans.copy(o);
	return ans;
}

Trie& Trie::operator = (const Trie& o) {
	cnt = o.cnt, tot_cnt = o.tot_cnt;
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		next[i] = NULL;
	}
	return *this;
}

Trie& Trie::operator += (const Trie& o) {
	this->copy(o);
	return *this;
}

void Trie::copy(const Trie& o) {
	this->cnt += o.cnt;
	this->tot_cnt += o.tot_cnt;
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		if (o.next[i]) {
			if (this->next[i]) {
				this->next[i]->copy(*o.next[i]);
			} else {
				this->next[i] = new Trie();
				this->next[i]->copy(*o.next[i]);
			}
		}
	}
}

int Trie::get_cnt() const {
	return cnt;
}

int Trie::get_tot_cnt() const {
	return tot_cnt;
}

char* Trie::insert(char* s) {
	char c = *s;
	this->tot_cnt++;
	if (c >= 'A' && c <= 'Z') {
		c = c - 'A' + 'a';
	} else if (c < 'a' || c > 'z') {
		this->cnt++;
		return s;
	}
	int x = c - 'a';
	if (this->next[x]) {
		return this->next[x]->insert(s + 1);
	} else {
		this->next[x] = new Trie();
		return this->next[x]->insert(s + 1);
	}
}

int Trie::query(char* s) const {
	char c = *s;
	if (c >= 'A' && c <= 'Z') {
		c = c - 'A' + 'a';
	} else if (c < 'a' || c > 'z') {
		return this->get_cnt();
	}
	int x = c - 'a';
	if (this->next[x]) {
		return this->next[x]->query(s + 1);
	} else {
		return 0;
	}
}

int Trie::query_pre(char* s) const {
	char c = *s;
	if (c >= 'A' && c <= 'Z') {
		c = c - 'A' + 'a';
	} else if (c < 'a' || c > 'z') {
		return this->get_tot_cnt();
	}
	int x = c - 'a';
	if (this->next[x]) {
		return this->next[x]->query(s + 1);
	} else {
		return 0;
	}
}

int Trie::remove(char* s) {
	char c = *s;
	if (c >= 'A' && c <= 'Z') {
		c = c - 'A' + 'a';
	} else if (c < 'a' || c > 'z') {
		int res = this->cnt;
		this->tot_cnt -= res;
		this->cnt = 0;
		return res;
	}
	int x = c - 'a';
	if (this->next[x]) {
		int res = this->remove(s + 1);
		this->tot_cnt -= res;
		return res;
	} else {
		return 0;
	}
}

vector< pair<int, string> > Trie::get_all() const {
	vector< pair<int, string> > ans;
	ans.clear();
	get_dfs(0, ans);
	return ans;
}

void Trie::get_dfs(int deep, vector< pair<int, string> >& ans) const {
	static char s[1024];
	if (deep >= 1024) {
		return;
	}
	if (this->cnt) {
		s[deep] = '\0';
		string tmp = s;
		ans.push_back(make_pair(this->cnt, tmp));
	}
	for (int i = 0; i < CHAR_SET_NUM; i++) {
		if (this->next[i]) {
			s[deep] = 'a' + i;
			this->next[i]->get_dfs(deep + 1, ans);
		}
	}
}

