//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#include <algorithm>
#include "itensor/util/input.h"
#include "itensor/util/print.h"

using std::ostream;
using std::istream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

namespace itensor {


void InputFile::open()
    {
    close();
    file_.open(filename_.c_str());
    if(!file_) 
        error("can't open input file");
    opened_ = true;
    }

void InputFile::close()
    {
    if(opened_) 
        { 
        file_.clear(); 
        file_.close(); 
        }
    opened_ = false;
    }

ostream& 
operator<<(ostream &s, InputFile const& a)
    {
    auto fname = a.filename();
    InputFile f(fname);
    f.open();
    s << "Input filename is " << f.filename() << endl;
    char c;
    while(f.file().get(c)) s << c;
    f.close();
    return s;
    }

void eatwhite(istream& is)
    {
    char c;
    while(is.get(c))
	{
	if(!isspace(c))
	    {
	    is.putback(c);
	    break;
	    }
	}
    }

void matchbracket(istream& is)
    {
    int nest = 1;
    char c;
    while(is >> c)
	{
	if(c == '{')
	    nest++;
	else if(c == '}')
	    nest--;
	// cout << "match: " << c << endl;
	if(nest == 0) break;
	}
    if(nest != 0)
	error("unterminated bracket");
    }

int gettoken(istream& is, string& s)
    {
    char c[512], ch;
    eatwhite(is);
    for(int i = 0; i < 512; i++)
	{
	if(!is.get(ch))
	    return 0;
	// cout << "gettoken: " << ch << endl;
	if(isalpha(ch) || ch == '_' || (i != 0 && isdigit(ch)))
	    c[i] = ch;
	else
	    {
	    if(i == 0)
		{
		c[0] = ch;
		c[1] = '\0';
		s = c;
		return 1;
		}
	    is.putback(ch);
	    c[i] = '\0';
	    break;
	    }
	}
    s = c;
    return 1;
    }

InputGroup::
InputGroup(std::string filename, 
           std::string groupname,
           const char* c)
    : parent(0), 
      name_(groupname), 
      quiet(false)
    {
    infile.set_managed(new InputFile(filename));
    //std::cout << "Making input group " << name_;
    //if(c) std::cout << ": " << c;
    //std::cout << std::endl;
    }

InputGroup::
InputGroup(InputFile& inf, 
           std::string name,
           const char* c)
    : parent(0), 
      name_(name), 
      quiet(false)
    {
    infile.set_external(&inf);
    //std::cout << "Making input group " << name_;
    //if(c) std::cout << ": " << c;
    //std::cout << std::endl;
    }

InputGroup::
InputGroup(InputGroup& par, 
           std::string name,
           const char* c)
    : parent(&par), 
      name_(name),
      quiet(false)
    {
    infile.set_external(&(*par.infile));
    //std::cout << "Making input group " << parent->name_ << "." << name_;
    //if(c) std::cout << ": " << c;
    //std::cout << std::endl;
    }

InputGroup::
~InputGroup()
    {
    }

int InputGroup::GotoGroup()
    {
    if(parent != 0)
	parent->GotoGroup();
    else
	infile->open();
    eatwhite(infile->file());
    while(1)
	{
	string s;
	if(!gettoken(infile->file(),s))
	    return 0;
	// cout << "GotoGroup: got string " << s << endl;
	if(s == name_)
	    {
	    eatwhite(infile->file());
	    char c;
	    infile->file().get(c);
	    if(c != '{')
		error("bracket does not follow group name");
	    eatwhite(infile->file());
	    return 1;
	    }
	else
	    {
	    if(s[0] == '{')
		{
		matchbracket(infile->file());
		}

	    while(1)
		{
		char c;
		infile->file().get(c);
		if(c == '{')
		    matchbracket(infile->file());
		if(c == '\n')
		    {
		    eatwhite(infile->file());
		    break;
		    }
		}
	    }
	}
    }

int InputGroup::GotoToken(string s)
    {
    if(!GotoGroup()) 
	return 0;
    while(1)
	{
	string t;
	if(!gettoken(infile->file(),t))
	    {
	    // cout << "no token, returning from GotoToken" << endl;
	    return 0;
	    }
	// cout << "GotoToken: got string " << t << endl;
	if(t[0] == '}') return 0;
	if(t == s) 
	    {
	    eatwhite(infile->file());
	    char c;
	    if(!(infile->file().get(c))) return 0;
	    if(c != '=')
		return 0;
	    eatwhite(infile->file());
	    return 1;
	    }
	if(t[0] == '{')
	    {
	    matchbracket(infile->file());
	    }
	while(1)
	    {
	    char c;
	    if(!(infile->file().get(c))) return 0;
	    if(c == '\n' || c == ',' || c == ';')
		{
		eatwhite(infile->file());
		break;
		}
	    if(c == '{')
		matchbracket(infile->file());
	    }
	}
    }

int InputGroup::
GetInt(string s, 
       int& res,
       bool hasdf,
       int df)
    {
    if(!GotoToken(s) || !(infile->file() >> res)) 
        {
        if(!quiet && hasdf) printfln("Def %s.%s = %s",name_,s,df);
        return 0;
        }
    if(!quiet) printfln("Got %s.%s = %s",name_,s,res);
    return 1;
    }

int InputGroup::
GetLong(string s, 
        long& res,
        bool hasdf,
        long df)
    {
    if(!GotoToken(s) || !(infile->file() >> res)) 
        {
        if(!quiet && hasdf) printfln("Def %s.%s = %s",name_,s,df);
        return 0;
        }
    if(!quiet) printfln("Got %s.%s = %s",name_,s,res);
    return 1;
    }

int InputGroup::
GetReal(string s, 
        Real& res,
        bool hasdf,
        Real df)
    {
    if(!GotoToken(s) || !(infile->file() >> res)) 
        {
        if(!quiet && hasdf) printfln("Def %s.%s = %s",name_,s,df);
        return 0;
        }
    if(!quiet) printfln("Got %s.%s = %s",name_,s,res);
    return 1;
    }

int InputGroup::
GetString(string s, 
          string& res,
          bool hasdf,
          string df)
    {
    if(!GotoToken(s) || !(infile->file() >> res)) 
        {
        if(!quiet && hasdf) printfln("Def %s.%s = %s",name_,s,df);
        return 0;
        }
    if(!quiet) printfln("Got %s.%s = %s",name_,s,res);
    return 1;
    }

char mydolower(char c) { return tolower(c); }

int InputGroup::
GetYesNo(string s, 
         int& yes,
         bool hasdf,
         int df)
    {
    string res;
    if(!GotoToken(s) || !(infile->file() >> res)) 
        {
        if(!quiet && hasdf) printfln("Def %s.%s = %s",name_,s,df);
        return 0;
        }
    if(!quiet) printfln("Got %s.%s = %s",name_,s,res);
    transform(res.begin(),res.end(),res.begin(),mydolower);
    if(res == "yes" || res == "y" || res == "true")
        {
        yes = 1;
        return 1;
        }
    if(res == "no" || res == "n" || res == "false")
        {
        yes = 0;
        return 1;
        }
    return 0;
    }

int InputGroup::
GetYesNo(string s, 
         bool& yes,
         bool hasdf,
         bool df)
    {
    int resi = 0;
    int got = GetYesNo(s,resi,hasdf,df);
    yes = (resi==1);
    return got;
    }

void InputGroup::SkipLine()
    {
    char c = '\0';
    while(c != '\n')
	infile->file().get(c);
    eatwhite(infile->file());
    }

int InputGroup::
getInt(std::string s, int def)
    {
    int res = 0;
    int got = GetInt(s,res,true,def);
    if(!got) return def;
    return res;
    }

Real InputGroup::
getReal(std::string s, Real def)
    {
    Real res = 0;
    int got = GetReal(s,res,true,def);
    if(!got) return def;
    return res;
    }

std::string  InputGroup::
getString(std::string s, std::string def)
    {
    std::string res;
    int got = GetString(s,res,true,def);
    if(!got) return def;
    return res;
    }

bool  InputGroup::
getYesNo(std::string s, bool def)
    {
    bool res = false;
    int got = GetYesNo(s,res,true,def);
    if(!got) return def;
    return res;
    }

int InputGroup::
getInt(std::string s)
    {
    int res = 0;
    GetIntM(s,res);
    return res;
    }

Real InputGroup::
getReal(std::string s)
    {
    Real res = 0;
    GetRealM(s,res);
    return res;
    }

std::string  InputGroup::
getString(std::string s)
    {
    std::string res;
    GetStringM(s,res);
    return res;
    }

bool  InputGroup::
getYesNo(std::string s)
    {
    int res = false;
    GetYesNoM(s,res);
    return bool(res);
    }

void InputGroup::
GetIntM(string s, 
        int& res)
    {
    if(!GetInt(s,res))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::
GetLongM(string s, long& res)
    {
    if(!GetLong(s,res))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::
GetRealM(string s, Real& res)
    {
    if(!GetReal(s,res))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::
GetStringM(string s, string& res)
    {
    if(!GetString(s,res))
        error("mandatory item: " + s + ", exiting");
    }

void InputGroup::
GetYesNoM(string s, int& yes)
    {
    if(!GetYesNo(s,yes))
        error("mandatory item: " + s + ", exiting");
    }

} //namespace itensor

