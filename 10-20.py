
# coding: utf-8

# とりあえず放線菌属について属をたくさん抜き出してみる。

# In[5]:

import lxml.html
import requests


# In[3]:

html = requests.get("https://ja.wikipedia.org/wiki/放線菌")


# In[4]:

print(html)


# In[6]:

dom = lxml.html.fromstring(html.text)


# In[9]:

print(dom.xpath('//*[@id="mw-content-text"]/div/ul[9]/li/ul/li[1]/i[2]'))


# 一番下のiタグが科名+属名  
# 二番下のliタグが目内での順番  
# 三番下のulタグは目名とそれ以外の判別用  
# 四番下のliタグは綱内での目数  
# 五番下のulはwikiの上から順番

# In[10]:

print(dom)


# In[22]:

print(dom.xpath('//*[@id="mw-content-text"]/div/ul[9]')[0].xpath('li/ul/li[1]/i[2]')[0].text)


# In[23]:

(/*[@id="mw-content-text"]/div/ul[14])


# 9から14まで綱があり、それらの子要素をどんどん抜き出していく。

# In[25]:

for i in range(9, 15):
    print(i)


# In[50]:

for i1 in range(9, 15):
    moku = dom.xpath('//*[@id="mw-content-text"]/div/ul[' + str(i1) + ']/li/ul/li')
    # print(list(moku[0]))
    for i2 in moku:
        for i3 in list(i2)[1:]:
            if i3.tag == "a":
                print(list(i3)[0].text)
            elif i3.tag == "i":
                print(i3.text)
            else:
                print("DAME")
    break


# In[51]:

acti = []
for i1 in range(9, 15):
    moku = dom.xpath('//*[@id="mw-content-text"]/div/ul[' + str(i1) + ']/li/ul/li')
    # print(list(moku[0]))
    for i2 in moku:
        for i3 in list(i2)[1:]:
            if i3.tag == "a":
                acti.append(list(i3)[0].text)
            elif i3.tag == "i":
                acti.append(i3.text)
            else:
                print("DAME")


# In[52]:

print(acti)


# In[53]:

print(len(acti), len(set(acti)))


# 放線菌の属の抜き出しは終わった。  
# この程度の全抜きなら30分くらいで終わる。

# In[63]:

class MCS_Finder2(object):
    def __init__(self, genus):
        import urllib.request
        import os
        import time
        from rdkit.Chem import rdFMCS
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem import rdDepictor
        import requests
        import lxml.html
        import re
        import networkx as nx

        self.genus = genus
        self.mol_list = None

        if not os.path.exists("test/" + self.genus):
            os.mkdir("test/" + self.genus)

    def get_html(self, genus):
        """
        get htmlfile from KNApSAck search engine

        input
            genus: str, genusric name

        output
            html: requests.models.Response
        """
        import requests
        html = requests.get("http://kanaya.naist.jp/knapsack_jsp/result.jsp?sname=organism&word=" + genus)
        return html

    def get_Cnumber(self, html, limit=2000):
        """
        get Cnumber list from KNApSAck htmlfile

        input
            html: requests.models.Response
            limit: int, itertion limit

        output
            Cnumber: list, sorted list of Cnumber
        """
        import lxml.html
        dom = lxml.html.fromstring(html.text)
        i = 1
        Cnumber = set()
        genus = dom.xpath('//*[@id="my_contents"]/font[2]')[0].text
        genus = genus[0].upper() + genus[1:]
        while(True):
            if i > limit:
                print("max itertion")
                raise Exception("max iteration change limit")

            try:
                Cn = dom.xpath('//*[@class="sortable d1"]/tr[' + str(i) + ']/td[1]/a')[0].text
            except IndexError:
                # print("finish getting Cnumber")
                # print(i)
                break

            try:
                if genus != dom.xpath('//*[@class="sortable d1"]/tr[' + str(i) + ']/td[6]/font')[0].text:
                    i += 1
                    continue
            except IndexError:
                print("font error line ", + str(i))
                i += 1
                continue

            Cnumber.add(Cn)
            i += 1
        Cnumber = list(sorted(Cnumber))
        return Cnumber

    def search(self, Cnumberlist, filename=None):
        "make kcffile only given Cnumbers from all kcffile"
        if filename is None:
            filename = self.genus + "/kcfs.kcfs"
        with open(filename, "w"):
            pass

        for Cnumber in Cnumberlist:
            Cn = int(Cnumber[1:])
            if Cn <= 9217:
                page = "1-9"
            elif Cn <= 19275:
                page = "10-19"
            elif Cn <= 29326:
                page = "20-29"
            elif Cn <= 39355:
                page = "30-39"
            elif Cn <= 49370:
                page = "40-49"
            elif Cn <= 50409:
                page = "50-59"
            else:
                print("your Cnumber " + Cnumber + " is too large")
                continue

            with open("../../../database/kcfs/KNApSAck" + page + ".kcfs") as f:
                for (i, line) in enumerate(f):
                    # print(line[12:21]) # C00000000
                    if line[12:21] == Cnumber:
                        # temp = i
                        # print("find", temp)
                        lin = line
                        flag = 0
                        with open(filename, "a") as fw:
                            while((lin[:5] != "ENTRY" and lin != "") or flag != 1):
                                flag = 1
                                fw.write(lin)
                                lin = f.readline()
                            else:
                                break
        return True

    def kcfs2count(self, kcfs, txt):
        import re
        dic = {}
        with open(kcfs, "r") as f:
            for mol in f.read(None).split("///\n"):
                # sta = 0
                sta2 = 0
                type_ = 0
                for line in mol.split("\n"):
                    # if re.match("^\S", line):
                        # sta = line.split()[0]
                    if re.match("^\s\s\S", line):
                        sta2 = line.split()[0]
                    if re.match("///", line):
                        pass
                    elif sta2:
                        a = line[12:].split()
                        type_ = sta2
                        try:
                            num = int(re.findall("\d+", a[1])[0])
                            str_ = a[0]
                        except IndexError:
                            continue
                        str1 = re.sub("[a-z]", "", str_)
                        str2 = re.sub("\d", "", str1)
                        dic.setdefault((type_, str_), 0)
                        dic[(type_, str_)] += num
                        dic.setdefault((type_, str2), 0)
                        dic[(type_, str2)] += num
                        dic.setdefault((type_, str1), 0)
                        dic[(type_, str1)] += num

        array = []
        for item in dic.items():
            array += [[0 - item[1], item[0]]]

        with open(txt, "w") as f2:
            index = 0
            for list_ in sorted(array):
                index += 1
                num = str(index)
                while(len(num) < 8):
                    num = "0" + num
                num = "S" + num + list_[1][0][0]
                f2.write(num + "\t" + list_[1][0] + "\t" + list_[1][1] + "\t" + str(0 - list_[0]) + "\n")
        return True

    def split(self, countfile, result="test.txt", limit=0):
        with open(countfile, "r") as f:
            with open(result, "w") as f2:
                for line in f:
                    i = -1
                    if line[11] == "R":
                        i = -1
                        while(True):
                            if line[i] == "\t":
                                # print(line[i+1:])
                                break
                            i -= 1
                        if int(line[i + 1:]) < limit:
                            break
                        # print(line)
                        f2.write(line)
                    elif line[11] == "I":
                        i = -1
                        while(True):
                            if line[i] == "\t":
                                # print(line[i+1:])
                                break
                            i -= 1
                        if int(line[i + 1:]) < limit:
                            break
                        # print(line)
                        f2.write(line)
                    elif line[11] == "S":
                        i = -1
                        while(True):
                            if line[i] == "\t":
                                # print(line[i+1:])
                                break
                            i -= 1
                        if int(line[i + 1:]) < limit:
                            break
                        # print(line)
                        f2.write(line)
        return True


# In[66]:

def make_kcfs2(f, limit=2000, splimit=0):
    html = f.get_html(f.genus)
    Cnumber = f.get_Cnumber(html, limit)
    if Cnumber == []:
        return False
    f.search(Cnumber, "test/" + f.genus + "/kcfs.kcfs")
    f.kcfs2count("test/" + f.genus + "/kcfs.kcfs", "test/" + f.genus + "/kcfscount.txt")
    f.split("test/" + f.genus + "/kcfscount.txt", "test/" + f.genus + "/splitedcount.txt", splimit)
    return True


# In[78]:

get_ipython().run_cell_magic('time', '', 'import time\nimport datetime\nimport os\nprint(datetime.datetime.now())\nfor i in acti:\n    if os.path.exists("test/" + i):\n        continue\n    f = MCS_Finder2(i)\n    make_kcfs2(f)\n    time.sleep(10)')


# In[73]:

import datetime
print(datetime.datetime.now())


# In[ ]:

32

