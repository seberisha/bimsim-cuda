#ifndef RTS_ARGUMENTS
#define RTS_ARGUMENTS

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <algorithm>

#ifdef _WIN32
#include <Windows.h>
#endif

namespace rts{

class argument
{
private:
    bool ansi;

    //argument name
    std::string name;

    //description of the argument
    std::vector<std::string> desc;

    //argument values
    std::vector<std::string> vals;

    //range or example
    std::string range;

    //flag is true when the argument is user-specified
    bool flag;

    void parse_val(const std::string &s){

        vals.clear();

        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, ' ')) {
            vals.push_back(item);
        }
    }

    void parse_desc(const std::string &s){

        desc.clear();

        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, '\n')) {
            desc.push_back(item);
        }
    }

public:
    void set_ansi(bool b){ ansi = b; }
    //create an argument with a given name, description, and default value
    argument(std::string _name, std::string _desc, std::string _default = "", std::string _range = "")
    {
        name = _name;
        parse_desc(_desc);
        parse_val(_default);

        //if a default value is provided, set the flag
        if(_default != "")
            flag = true;
        else flag = false;

        range = _range;


    }

    int nargs()
    {
        return vals.size();
    }

    //return the value of a text argument
    std::string as_text(int n = 0)
    {
        if(!flag)
        {
            std::cout<<"ERROR - Argument requested without being set: "<<name<<std::endl;
            exit(1);
        }

        if(vals.size() > n)
            return vals[n];

        else return "";
    }

    //return the value of a floating point argument
    float as_float(int n = 0)
    {
        if(!flag)
        {
            std::cout<<"ERROR - Argument requested without being set: "<<name<<std::endl;
            exit(1);
        }

        if(vals.size() > n)
        {
            float r;
            if ( ! (istringstream(vals[n]) >> r) ) r = 0;
            return r;
        }

        else return 0;
    }

    //return the value of an integer argument
    int as_int(int n = 0)
    {
        if(!flag)
        {
            std::cout<<"ERROR - Argument requested without being set: "<<name<<std::endl;
            exit(1);
        }

        if(vals.size() > n)
        {
            int r;
            if ( ! (istringstream(vals[n]) >> r) ) r = 0;
            return r;
        }

        else return 0;
    }

    //get the width of the left column
    int col_width()
    {
        int n = 3;
        //add the length of the argument name
        n += name.size();

        //if there are any default parameters
        if(vals.size() > 0)
        {
            //padding (parenthesis, =, etc.)
            n += 6;

            //for each default argument value
            for(int v=0; v<vals.size(); v++)
                n += vals[v].size() + 1;
        }

        //add a buffer of 4 characters
        n += 4;

        return n;
    }


    //string output
    std::string toStr(int width = 0)
    {
        std::stringstream ss;

        int color_size = 0;


        //create the left column
        std::string left_part = std::string(" --") + name;
        if(vals.size() != 0)
        {
            if(ansi)
                left_part += "\033[1;32m";
            left_part += " ( = ";
            for(int v=0; v<vals.size(); v++)
                left_part += vals[v] + std::string(" ");
            left_part += ")";
            if(ansi)
                left_part += "\033[0m";
            if(ansi)
                color_size = 11;
        }
        else
            color_size = 0;

        //if no width is passed, put 4 spaces between left and right columns
        if(width == 0) width = col_width();

        ss<<std::left<<std::setw(width + color_size)<<left_part;

        //output right column
        for(int d=0; d<desc.size(); d++)
        {
            if(d == 0)
                ss<<desc[0];
            else
                ss<<std::endl<<setfill(' ')<<setw(width)<<" "<<desc[d];

        }

        //output the range in the right column
        if(range != "" && ansi)
            ss<<std::endl<<setfill(' ')<<setw(width)<<" "<<"    "<<std::string("\033[1;36m") + range + "\033[0m";
        else if(range != "")
            ss<<std::endl<<setfill(' ')<<setw(width)<<" "<<"    "<<range;

        return ss.str();
    }

    //compare the name of the argument to a string
    bool operator==(std::string rhs)
    {
        return (name == rhs);
    }

    //set the argument to a given value
    void set(std::string _value)
    {
        parse_val(_value);

        //set the flag
        flag = true;
    }

    bool is_set()
    {
        return flag;
    }

};

struct argsection
{
    std::string name;
    unsigned int index;
};

class arglist
{
private:
    bool ansi;

    //vector of arguments
    std::vector<argument> args;

    //column width of the longest argument
    int col_width;

    //list of sections
    std::vector<argsection> sections;

public:

    arglist(){ col_width = 0; }

    void set_ansi(bool b)
    {
        ansi = b;
        for(int i=0; i<args.size(); i++)
            args[i].set_ansi(ansi);
    }

    void add(std::string _name, std::string _desc, std::string _default = "", std::string _range = "")
    {
        argument arg(_name, _desc, _default, _range);
        arg.set_ansi(ansi);
        args.push_back(arg);

        col_width = std::max<int>(col_width, arg.col_width());
    }

    void section(std::string _name)
    {
        argsection s;
        s.name = _name;
        s.index = args.size();
        sections.push_back(s);
    }

    //output the arguments (generally in response to --help)
    std::string toStr()
    {
        std::stringstream ss;

        int si = -1;

        if(sections.size() > 0)
            si = 0;

        //for each argument
        for(int a=0; a<args.size(); a++)
        {
            if(si != -1 && a == sections[si].index)
            {
                if(ansi)
                    ss<<std::endl<<std::left<<setfill('=')<<setw(col_width)<<std::string("\033[1;31m") + sections[si].name<<"\033[0m"<<std::endl;
                else
                    ss<<std::endl<<std::left<<setfill('=')<<setw(col_width)<<sections[si].name<<std::endl;
                si++;
                if(si == sections.size()) si = -1;
            }

            ss<<args[a].toStr(col_width)<<std::endl;
        }

        return ss.str();
    }

    int index(std::string _name)
    {
        int i = find(args.begin(), args.end(), _name) - args.begin();

        if(i >= args.size())
            i = -1;

        return i;
    }

    void set(std::string _name, std::string _value)
    {
        int i = index(_name);

        if(i != -1)
        {
            args[i].set(_value);
            //adjust the column width if necessary
            col_width = max(col_width, args[i].col_width());
        }
        else
            std::cout<<"ERROR - Argument not recognized: "<<_name<<std::endl;
    }

    //parse a parameter string
    void parse(int argc, char* argv[])
    {
        //if the number of arguments is 1, we're done
        if(argc <= 1) return;

        std::string name;
        std::string params;

        for(int i=1; i<argc; i++)
        {
            //if the argument is a parameter name
            if(argv[i][0] == '-' && argv[i][1] == '-')
            {
                //add any previous arguments
                if(name != "")
                    set(name, params);
                //set the current argument to this name
                name = argv[i]+2;

                //std::cout<<"\t\t name: "<<name<<std::endl;

                //clear the parameters list
                params = "";
            }
            else
            {
                if(params != "")
                    params += " ";
                params += argv[i];
            }
        }

        //set the last argument
        set(name, params);
    }

    //determine if a parameter has been set (either specified by the user or with a default value)
    bool operator()(std::string _name)
    {
       // std::cout<<"\t\t argslist: _name: "<<_name<<std::endl;
        int i = find(args.begin(), args.end(), _name) - args.begin();

        if(i < 0)
        {
            std::cout<<"ERROR - Unspecified parameter name: "<<_name<<std::endl;
            exit(1);
        }

        return args[i].is_set();
    }

    int nargs(std::string _name)
    {
        int i = find(args.begin(), args.end(), _name) - args.begin();

        if(i < 0)
        {
            std::cout<<"ERROR - Unspecified parameter name: "<<_name<<std::endl;
            exit(1);
        }

        return args[i].nargs();
    }

    argument operator[](std::string _name)
    {
        int i = find(args.begin(), args.end(), _name) - args.begin();

        if(i < 0)
        {
            std::cout<<"ERROR - Unspecified parameter name: "<<_name<<std::endl;
            exit(1);
        }
        return args[i];
    }


};



}	//end namespace rts

#endif
