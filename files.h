#ifndef _SHANGKE_FILES_H_
#define _SHANGKE_FILES_H_

#include <vector>
#include <string>
#include <string.h>
#include <io.h>

using namespace std;

class Files: public vector<string> {
	public:
		/* *
		 * 得到一个目录下的所有文件
		 */
		void getFiles(string path, bool needclear = true) {
			if (needclear) {
				this->clear();
			}
			int hFile = 0;
			struct _finddata_t fileinfo;
			string p;
			if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1) {
				do {
					if (fileinfo.attrib & _A_SUBDIR) {
						if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0) {
							getFiles(p.assign(path).append("\\").append(fileinfo.name), false);
						}
					} else {
						push_back(p.assign(path).append("\\").append(fileinfo.name));
					}
				} while (_findnext(hFile, &fileinfo) == 0);
				_findclose(hFile);
			}
		}
};

#endif
