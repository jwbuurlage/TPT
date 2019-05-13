#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace tpt {
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
        ss << result;
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

    void print(std::ostream& os = std::cout) {
        os << title_ << "\n";

        std::string hline = "|";
        auto repeat = [](std::string base, auto n) -> std::string {
            std::string result = "";
            while (n--) {
                result += base;
            }
            return result;
        };
        hline += repeat("-", row_size_ + 2);
        hline += "|";
        for (auto col : columns_) {
            hline += repeat("-", column_width_[col] + 2);
            hline += "|";
        }

        auto addElement = [repeat](int width, std::stringstream& result,
                             std::string entry) {
            result << entry;
            result << repeat(" ", width - entry.size());
        };

        std::stringstream ss;
        ss << "| ";
        addElement(row_size_, ss, row_title_);

        for (auto& col : columns_) {
            ss << " | ";
            addElement(column_width_[col], ss, col);
        }
        ss << " |";

        os << hline << "\n";
        os << ss.str() << "\n";
        os << hline << "\n";

        for (auto& rowCols : entries_) {
            std::stringstream rowSs;
            rowSs << "| ";
            addElement(row_size_, rowSs, rowCols.first);
            for (auto& col : columns_) {
                rowSs << " | ";
                addElement(column_width_[col], rowSs, rowCols.second[col]);
            }
            rowSs << " |";
            os << rowSs.str() << "\n";
        }
        os << hline << "\n";
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
} // namespace tpt
