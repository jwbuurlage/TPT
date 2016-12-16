#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <string>

namespace tomo {
namespace util {

class report {
   public:
    report(std::string title, std::string rowTitle)
        : title_(title), rowTitle_(rowTitle) {
        rowSize_ = rowTitle_.size();
    }

    void add_column(std::string colName, std::string texName = "") {
        columns_.push_back(colName);
        if (!texName.empty())
            columnsTex_.push_back(texName);
        else
            columnsTex_.push_back(colName);
        columnWidth_[colName] = colName.size();
    }

    void add_row(std::string row) {
        entries_[row] = std::map<std::string, std::string>();

        if (row.size() > rowSize_) {
            rowSize_ = row.size();
        }
    }

    template <typename T>
    void add_result(std::string row, std::string column, T result) {
        if (entries_.find(row) == entries_.end()) {
            std::cout << "Error: Trying to add result to non-existing row\n";
            return;
        }
        std::stringstream ss;
        ss << std::fixed << std::setprecision(1) << result;
        entries_[row][column] = ss.str();
        entriesTex_[row][column] = ss.str();

        if (ss.str().size() > columnWidth_[column]) {
            columnWidth_[column] = ss.str().size();
        }
    }

    void add_result(std::string row, std::string column, std::string result,
                   std::string texResult) {
        add_result(row, column, result);
        entriesTex_[row][column] = texResult;
    }

    void print() {
        std::cout << title_ << "\n";

        unsigned int lineSize = rowSize_ + 4;
        for (auto col : columnWidth_) {
            lineSize += col.second + 2;
        }
        std::string hline = "";
        for (unsigned int i = 0; i < lineSize; ++i) hline.push_back('-');

        auto addElement = [](int width, std::stringstream& result,
                             std::string entry) {
            result << std::left << std::setprecision(1) << std::setw(width)
                   << std::setfill(' ') << entry;
        };

        std::stringstream ss;
        addElement(rowSize_ + 2, ss, rowTitle_);
        ss << "| ";

        for (auto& col : columns_) addElement(columnWidth_[col] + 2, ss, col);

        std::cout << hline << "\n";
        std::cout << ss.str() << "\n";
        std::cout << hline << "\n";

        for (auto& rowCols : entries_) {
            std::stringstream rowSs;
            addElement(rowSize_ + 2, rowSs, rowCols.first);
            rowSs << "| ";
            for (auto& col : columns_) {
                addElement(columnWidth_[col] + 2, rowSs, rowCols.second[col]);
            }
            std::cout << rowSs.str() << "\n";
        }
        std::cout << hline << "\n";
    }

    void save_to_csv();
    void read_from_csv();

    void save_to_tex(std::string filename) {
        auto replaceTex = [](std::string entry) {
            std::string texEntry = entry;

            auto pos = entry.find("+-");
            if (pos != std::string::npos) {
                texEntry =
                    texEntry.substr(0, pos) + "\\pm" + texEntry.substr(pos + 2);
            }

            pos = entry.find("%");
            if (pos != std::string::npos) {
                texEntry =
                    texEntry.substr(0, pos) + "\\%" + texEntry.substr(pos + 1);
            }

            return texEntry;
        };

        std::ofstream fout(filename);
        fout << "\\begin{table}" << std::endl;
        fout << "\\centering" << std::endl;
        fout << "\\begin{tabular}{|l|";
        for (unsigned int i = 0; i < columns_.size(); ++i) {
            fout << "l";
            fout << ((i < (columns_.size() - 1)) ? " " : "|}");
        }
        fout << std::endl << "\\hline" << std::endl;
        fout << "\\textbf{" << rowTitle_ << "} & ";
        for (unsigned int i = 0; i < columns_.size(); ++i) {
            fout << "$" << columnsTex_[i] << "$";
            fout << ((i < (columns_.size() - 1)) ? " & " : "\\\\");
        }
        fout << std::endl;
        fout << "\\hline" << std::endl;

        for (auto& rowCols : entriesTex_) {
            fout << "\\verb|" << rowCols.first << "| & ";
            for (unsigned int i = 0; i < columns_.size(); ++i) {
                fout << "$" << replaceTex(rowCols.second[columns_[i]]) << "$";
                fout << ((i < (columns_.size() - 1)) ? " & " : "\\\\");
            }
            fout << std::endl;
        }
        fout << "\\hline" << std::endl;
        fout << "\\end{tabular}" << std::endl;
        ;

        fout << "\\caption{\\ldots}" << std::endl;
        fout << "\\end{table}" << std::endl;
    }

   private:
    std::string title_;
    std::string rowTitle_;
    std::map<std::string, std::map<std::string, std::string>> entries_;
    std::map<std::string, std::map<std::string, std::string>> entriesTex_;
    std::vector<std::string> columns_;
    std::vector<std::string> columnsTex_;
    std::map<std::string, unsigned int> columnWidth_;
    unsigned int rowSize_;
};

}  // namespace util
}  // namespace tomo
