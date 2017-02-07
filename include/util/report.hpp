#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace tomo {
namespace util {

class report {
  public:
    report(std::string title, std::string row_title)
        : title_(title), row_title_(row_title) {
        row_size_ = row_title_.size();
    }

    void add_column(std::string col_name, std::string tex_name = "") {
        columns_.push_back(col_name);
        if (!tex_name.empty())
            columns_tex_.push_back(tex_name);
        else
            columns_tex_.push_back(col_name);
        column_width_[col_name] = col_name.size();
    }

    void add_row(std::string row) {
        entries_[row] = std::map<std::string, std::string>();

        if (row.size() > row_size_) {
            row_size_ = row.size();
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
        entries_tex_[row][column] = ss.str();

        if (ss.str().size() > column_width_[column]) {
            column_width_[column] = ss.str().size();
        }
    }

    void add_result(std::string row, std::string column, std::string result,
                    std::string texResult) {
        add_result(row, column, result);
        entries_tex_[row][column] = texResult;
    }

    void print() {
        std::cout << title_ << "\n";

        unsigned int line_size = row_size_ + 4;
        for (auto col : column_width_) {
            line_size += col.second + 2;
        }
        std::string hline = "";
        for (unsigned int i = 0; i < line_size; ++i)
            hline.push_back('-');

        auto addElement = [](int width, std::stringstream& result,
                             std::string entry) {
            result << std::left << std::setprecision(1) << std::setw(width)
                   << std::setfill(' ') << entry;
        };

        std::stringstream ss;
        addElement(row_size_ + 2, ss, row_title_);
        ss << "| ";

        for (auto& col : columns_)
            addElement(column_width_[col] + 2, ss, col);

        std::cout << hline << "\n";
        std::cout << ss.str() << "\n";
        std::cout << hline << "\n";

        for (auto& rowCols : entries_) {
            std::stringstream rowSs;
            addElement(row_size_ + 2, rowSs, rowCols.first);
            rowSs << "| ";
            for (auto& col : columns_) {
                addElement(column_width_[col] + 2, rowSs, rowCols.second[col]);
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
        fout << "\\textbf{" << row_title_ << "} & ";
        for (unsigned int i = 0; i < columns_.size(); ++i) {
            fout << "$" << columns_tex_[i] << "$";
            fout << ((i < (columns_.size() - 1)) ? " & " : "\\\\");
        }
        fout << std::endl;
        fout << "\\hline" << std::endl;

        for (auto& rowCols : entries_tex_) {
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
    std::string row_title_;
    std::map<std::string, std::map<std::string, std::string>> entries_;
    std::map<std::string, std::map<std::string, std::string>> entries_tex_;
    std::vector<std::string> columns_;
    std::vector<std::string> columns_tex_;
    std::map<std::string, unsigned int> column_width_;
    unsigned int row_size_;
};

} // namespace util
} // namespace tomo
