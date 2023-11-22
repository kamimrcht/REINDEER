/******************************************************************************\
*                                                                              *
*  Copyright (C) 2023-      Anthony Boureux                                    *
*  Copyright (C) 2023       Nathan Boureux                                     *
*                                                                              *
*                                                                              *
*  This file is part of Reindeer program.                                      *
*                                                                              *
*  The purpose of Reindeer_server is to run Reindeer:                          *
*    - load an index                                                           *
*    - remains in memory until it receive a stop order                         *
*    - answer to query via a socket                                            *
*                                                                              *
*   Reinder_server is free software: you can redistribute it and/or modify     *
*   it under the terms of the GNU General Public License as published by       *
*   the Free Software Foundation, either version 3 of the License, or          *
*   (at your option) any later version.                                        *
*                                                                              *
*   Reinder_server is distributed in the hope that it will be useful,          *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
*   GNU General Public License for more details.                               *
*                                                                              *
*   You should have received a copy of the GNU General Public License          *
*   along with Reindeer program. If not, see <http://www.gnu.org/licenses/>.  *
*                                                                              *
*                                                                              *
*******************************************************************************/

#include <filesystem>
#include <sys/socket.h> // For socket functions
#include <netinet/in.h> // For sockaddr_in
#include <getopt.h>
#include <unistd.h>
#include "reindeer.hpp"
#include "utils.hpp"
#include "../version.h" // manage version from git

using namespace std;
namespace fs = filesystem;

typedef struct data_for_socket{
    int port = 0;
    string index_directory = "";
    string command_file = "";
    bool verbose = 0;
} DataSoc;

void help() {
    cout << endl;
    cout <<
    "* Mandatory parameters *\n"
    "  -p, --port <number>                :       Port\n"
    "* Client mode parameter *\n"
    "  -i, --input <file>                 :       Input file with command to send to a reindeer_socket server\n"
    "* Server mode parameter *\n"
    "  -l <directory>                     :       Reindeer index directory (should be reindeer_index_files if you've not used -o during indexing)\n"
    "\n"
    "* General parameters *\n"
    "  --verbose, -v                      :       Verbose mode\n"
    "  --version, -V                      :       Show version\n"
    "  -h, --help                         :       Print help (what you are currently seeing)\n"
    << endl;
    exit(0);
}

void process_args(int argc, char** argv, DataSoc& data) {

    const char* const short_opts = "vVhp:l:";
    const option long_opts[] = {
        { "port", required_argument, nullptr, 'p' },
        { "input", required_argument, nullptr, 'i' },
        { "index", required_argument, nullptr, 'l' },
        { "verbose", no_argument, nullptr, 'v' },
        { "version", no_argument, nullptr, 'V' },
        { "help", no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0 }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt) {
            break;
        }
        switch (opt) {
            case 'p':
                data.port = stoi(optarg);
                break;
            case 'i':
                data.command_file = optarg;
                break;
            case 'l':
                data.index_directory = optarg;
                break;
            case 'V':
                cout << VERSION << endl;
                exit(0);
                break;
            case 'v':
                data.verbose = 1;
                break;
            case 'h':
            case '?': // Unrecognized option
            default:
                help();
                break;
        }
    }
}

int main (int argc, char* argv[]) {

    if (argc == 1 || argc == 3 || argc == 4) {
        cout << strerror(EINVAL) << " : Missing argument(s) | Try -h or --help" << endl;
        exit(EINVAL);
    }

    DataSoc data;
    process_args(argc, argv, data);

    string message;

    // Create a socket (IPv4, TCP)
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1) {
        cerr << "Failed to create socket. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    if (data.verbose)
        cerr << " * Socket created\n";

    // Listen to port data.port on any address
    sockaddr_in sockaddr;
    sockaddr.sin_family = AF_INET;
    sockaddr.sin_addr.s_addr = INADDR_ANY;
    sockaddr.sin_port = htons(data.port); // htons is necessary to convert a number to

      // ********************
      // in client mode
      // ********************
    if (data.index_directory == "") {
      if (connect(sockfd, (struct sockaddr*) &sockaddr, sizeof(sockaddr)) < 0) {
        cerr << "Connection failed to server on port " << data.port << endl;
        exit(EXIT_FAILURE);
      }
      if (data.verbose)
        cerr << " * Connection established with the server on port " << data.port << endl;

      // wait welcome message
      auto bytesRead = read(sockfd, buffer, 256);
      if (strlen(buffer))
        cerr << " <- " << buffer << endl ;
      // ask for index
      message = "INDEX\n";
      write(sockfd, message.c_str(), message.size());
    }
    else {
      // ********************
      // in server mode
      // ********************

    // avoid err 98: socket already in use, due to timeout when server close
    int yes = 1;
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (void*)&yes, sizeof(yes)) < 0) {
        cerr << "setsockopt() failed. Error: errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    if (data.verbose)
        cerr << " * Connection options set to SOL_SOCKET and SO_REUSEADDR\n";
    // network byte order
    if (bind(sockfd, (struct sockaddr*)&sockaddr, sizeof(sockaddr)) < 0) {
        cerr << "Failed to bind to port " << data.port << ". errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    if (data.verbose)
        cerr << " * Connection bind to port " << data.port << endl;

    // Start listening. Hold at most 10 connections in the queue
    if (listen(sockfd, 1) < 0) {
        cerr << "Failed to listen on socket. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    if (data.verbose) {
       cerr << " * Going to listen message from client" << endl;
       cerr << " * Waiting for client to connect" << endl;
    }

    // Grab a connection from the queue
    auto addrlen = sizeof(sockaddr);
    int connection = accept(sockfd, (struct sockaddr*)&sockaddr, (socklen_t*)&addrlen);
    if (connection < 0) {
        cerr << "Failed to grab connection. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    if (data.verbose)
      cerr << " * connected to one client" << endl;
    // give a welcome message
    message = "WELCOME to Reindeer socket server\n";
    write(connection, message.c_str(), message.size());

    string output = "reindeer_index_files", query_file = "", format = "raw";
    // create reindeer index object
    Reindeer_Index<uint16_t> reindeer_index(data.index_directory, output, 1, true);
    //load index
    reindeer_index.load_index();

    // give which index is used
    message = "INDEX:";
    message.append(reindeer_index.matrix_name);
    message.append("\n");
    write(connection, message.c_str(), message.size());
    if (data.verbose)
      cerr << message << endl;

    // Read from the connection
    bool quit = false;

    while (!quit) {
        char buffer[256];
        memset(buffer, 0, 255);
        auto bytesRead = read(connection, buffer, 256);

        // check if something on connection
        // if something but len=0 => client hangout
        if (strlen(buffer) > 2) {
            // remove the \n send
            buffer[bytesRead] = '\0';
            string outFile;

            // convert buffer to string obj
            string strBuffer(buffer, bytesRead);
            // remove end char if control ASCII
            strBuffer.erase(strBuffer.find_last_not_of(" \t\n\r\f\v") + 1);
            // analysis args coming from socket message
            switch(toupper(strBuffer.at(0))) {
                // QUIT message: Quit and Send a message to the connection
                case 'Q':
                    message = "See you soon !\n";
                    cerr << message << endl;
                    write(connection, message.c_str(), message.size());
                    quit = true;
                    break;
                // HELP message
                case 'H':
                    message = " * HELP = help message\n";
                    message.append(" * QUIT = quit message\n");
                    message.append(" * INDEX = ask index in use\n");
                    message.append(" * FILE:myfile.fasta[:THRESHOLD:value][:OUTFILE:myoutfile][:FORMAT:format]\n");
                    write(connection, message.c_str(), message.size());
                    break;
                // INDEX message: ask which index is in memory
                case 'I':
                    message = "INDEX:";
                    message.append(reindeer_index.matrix_name);
                    message.append("\n");
                    write(connection, message.c_str(), message.size());
                    break;
                case 'F':
                // informations are by pair
                { // to be able to initialize subcmds var, case are like goto
                    string query_file {};
                    vector<string> subcmds = split_utils(strBuffer, ':');
                    if (not (subcmds.size() % 2)) {
                        for (auto it = subcmds.begin(); it < subcmds.end(); it++) {
                            char option = toupper(string(*it).at(1));
                            switch(option) {
                                // get fasta file: FILE
                                case 'I':
                                    it++;
                                    query_file = *it;
                                    cerr << "FILE = " << query_file << endl;
                                    break;
                                // get threshold: THRESHOLD
                                case 'H':
                                    it++;
                                    reindeer_index.threshold = stoi(*it);
                                    cerr << "THRESHOLD = " << reindeer_index.threshold  << endl;
                                    break;
                                // get output file: OUTFILE
                                case 'U':
                                    it++;
                                    outFile = *it;
                                    cerr << "OUTFILE = " << outFile << endl;
                                    break;
                                case 'O':
                                    it++;
                                    reindeer_index.output_format = *it;
                                    cerr << "FORMAT = " << reindeer_index.output_format << endl;
                                    break;
                                default :
                                    message = "UNKNOWN command: ";
                                    message.append(1, option);
                                    message.append("\n");
                                    cerr << message << endl;
                                    write(connection, message.c_str(), message.size());
                            }
                        }
                        if (fs::exists(query_file)) {
                            if (!outFile.empty()) reindeer_index.output = outFile;
                            reindeer_index.querying(query_file, reindeer_index.threshold, reindeer_index.output_format);
                            // we have finished
                            message = "DONE\n";
                            write(connection, message.c_str(), message.length());
                        } else {
                            message = "ERROR:The entry is not a file or not given (FILE:path.fa is required)\n";
                            cerr << message << endl;
                            write(connection, message.c_str(), message.size());
                        }
                    } else {
                        message = "ERROR: no file given to command FILE, or not paired command\n";
                        write(connection, message.c_str(), message.size());
                    }  // if subcmds.size %2
                    break;
                }
                default :
                    message = "UNKNOWN command: ";
                    message.append(1, strBuffer.at(0));
                    message.append("\n");
                    write(connection, message.c_str(), message.size());
            }
        } // if len > 2
    } // while
    // Close the connections
    close(connection);
    } // end of server mode
    close(sockfd);
    return 0;
}
