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

// ************
// Message Command available
namespace Command {
  // require for switch
  enum Type {
    HELP,
    QUIT,
    STOP,
    INDEX,
    FILECMD,
    THRESHOLD,
    OUTFILE,
    FORMAT,
    max_types
  };
  // define array to convert string to enum
  constexpr std::array<std::string_view, max_types> cmdName {"HELP", "QUIT", "STOP", "INDEX", "FILE", "THRESHOLD", "OUTFILE", "FORMAT"};
  // asset both enum and array are identical
  static_assert(std::size(cmdName) == max_types);
};

constexpr Command::Type getCmd(string_view cmd) {
  // find the Id
  for (std::size_t index = 0; index < Command::cmdName.size(); index++)
    // use substr for FILECMD that is not a uniq word
    if (cmd.substr(0,2) == Command::cmdName[index].substr(0,2))
      return static_cast<Command::Type>(index);
  return Command::max_types;
};
// ************

// *****************************************
// help message
const std::string help_command = R"(
Help for command to send to server

 * QUIT : the server is stopped
 * STOP : the client is stopped, not the server
 * HELP : this help
 * INDEX = ask index in use
 * FILE:myfile.fasta[:THRESHOLD:value][:OUTFILE:myoutfile][:FORMAT:format]
)";

void help() {
    cout << endl;
    cout <<
    "Usage:\n"
    " - To start a socket server:\n"
    "      reindeer_socket -p 5000 -l myReindeerIndex_Path\n"
    " - To start a socket client to communicate with the previous server:\n"
    "      reindeer_socket -p 5000 -i my_query.txt\n"
    "\n\n"
    "* Mandatory parameters *\n"
    "  -p, --port <number>                :       Port\n"
    "  + Client mode parameter:\n"
    "    -i, --input <file>                 :       Input file with command to send to a reindeer_socket server\n"
    "  + Server mode parameter:\n"
    "    -l <directory>                     :       Reindeer index directory (should be reindeer_index_files if you've not used -o during indexing)\n"
    "\n"
    "* General parameters *\n"
    "  --command                          :       Show command available between client and server\n"
    "  --verbose, -v                      :       Verbose mode\n"
    "  --version, -V                      :       Show version\n"
    "  -h, --help                         :       Print help (what you are currently seeing)\n"
    << endl;
    exit(0);
}
// *****************************************

// *****************************************
// Manage arguments

void process_args(int argc, char** argv, DataSoc& data) {

    const char* const short_opts = "p:i:l:cvVh";
    const option long_opts[] = {
        { "port", required_argument, nullptr, 'p' },
        { "input", required_argument, nullptr, 'i' },
        { "index", required_argument, nullptr, 'l' },
        { "command", no_argument, nullptr, 'c' },
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
            case 'c':
                cout << help_command << endl;
                exit(0);
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
// *****************************************

// *****************************************
// Send message to the socket
bool send_message(int socket, string message, bool verbose = 0) {
  // append "\n"
  //message.append("\n");
  if (verbose)
    cerr << "-> " << message << endl;
  if (send(socket, message.c_str(), message.size(), 0) < 0) {
    cerr << "Error when sending message: " << message << endl;
    exit(EXIT_FAILURE);
  }
  return 1;
}

// *****************************************

// *****************************************
// Function to wait for message from socket

string wait_message(int socket, bool verbose = 0) {
  bool waitMessage = true;
  char buffer[256];
  memset(buffer, 0, 255);
  ssize_t bytesRead = 0;

  while (waitMessage) {
    bytesRead = recv(socket, buffer, 256, 0);
    // cut the buffer to length receive
    buffer[bytesRead] = '\0';
    if (bytesRead < 0) {
      cerr << "Error when receiving message\n";
      return "Error";
    } else {
      if (verbose)
        cerr << " <- " << buffer << endl ;
      waitMessage = false;
    }
  }
  // convert buffer to string obj
  string strBuffer(buffer, bytesRead);
  // remove end char if control ASCII
  strBuffer.erase(strBuffer.find_last_not_of(" \t\n\r\f\v") + 1);
  return strBuffer;
}
// *****************************************

int main (int argc, char* argv[]) {

    DataSoc data;
    process_args(argc, argv, data);

    // check args
    // we require the port at least
    if (!data.port) {
      cerr << "Error, you need to give at least the port number to connect" << endl;
      help();
      exit(EINVAL);
    }


    char buffer[256];
    memset(buffer, 0, 255);

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
      message = wait_message(sockfd, data.verbose);
      if (message.at(0) == 'W') {
        if (data.verbose) {
          cerr << "Receive welcome message from server\n";
        }
      }

      // ask for index
      message = "INDEX";
      //send_message(sockfd, message);
      //message = wait_message(sockfd, data.verbose);
      cerr << "Index is: " << message << endl;

      if (data.command_file == "USER") {
        // nothing given on cmd line, ask user
        bool waitCommand{true};
        while (waitCommand) {
          message = "";
          cout << "Give command to send to server:\n";
          cin >> message;
          if (message.at(0) == 'S')
            waitCommand = false;
          if (data.verbose)
            cerr << "Ask: " << message << endl;
          send_message(sockfd, message, data.verbose);
          message = wait_message(sockfd, data.verbose);
        }
      } else {
            std:ifstream fp;
    std::istream &filein = fs::exists(data.command_file)
      ? [&]() -> std::istream& {fp.open(data.command_file); return fp;}()
      : std::cin;
    //unique_ptr< istream > filein;         // filehandle for command file
    //if (fs::exists(data.command_file))
    //  filein = unique_ptr< istream >(new ifstream(data.command_file)); //unique_ptr< istream >(new zstr::ifstream(tags_file));
    //else {
    //  // try to use stdin
    //  filein& = std::cin; //(new zstr::istream(cin));
    //}

    for (message; getline(filein, message); ) {
      if (data.verbose)
        cerr << "Ask: " << message << endl;
      send_message(sockfd, message, data.verbose);
      message = wait_message(sockfd, data.verbose);
    }
    if (fp.is_open()) {
      fp.close();
    }

    // send a STOP connection, to free server
    send_message(sockfd, "STOP", data.verbose);
    }
    }
    else {
      // ********************
      // in server mode
      // ********************

      // avoid err 98: socket already in use, due to timeout when server close
      int yes = 1;
      //if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (void*)&yes, sizeof(yes)) < 0) {
      //  cerr << "setsockopt() failed. Error: errno: " << errno << endl;
      //  exit(EXIT_FAILURE);
      //}
      //if (data.verbose)
      //  cerr << " * Connection options set to SOL_SOCKET and SO_REUSEADDR\n";

      // network byte order
      if (bind(sockfd, (struct sockaddr*)&sockaddr, sizeof(sockaddr)) < 0) {
          cerr << "Failed to bind to port " << data.port << ". errno: " << errno << endl;
          exit(EXIT_FAILURE);
      }
      if (data.verbose)
          cerr << " * Connection bind to port " << data.port << endl;

      // Start listening. Hold at most 10 connections in the queue
      if (listen(sockfd, 10) < 0) {
          cerr << "Failed to listen on socket. errno: " << errno << endl;
          exit(EXIT_FAILURE);
      }
      if (data.verbose) {
         cerr << " * Going to listen message from client" << endl;
         cerr << " * Waiting for client to connect" << endl;
      }

      // Loading index before client is connected
      if (data.verbose)
        cerr << "LOADING index " << data.index_directory << endl;

      string output = "reindeer_index_files", query_file = "", format = "raw";
      // create reindeer index object
      Reindeer_Index<uint16_t> reindeer_index(data.index_directory, output, 1, true);
      //load index
      reindeer_index.load_index();
      // validate loading to stdout
      cout << "LOADED index\n";

    // start listen on multiple connection
    bool quit = false;

    while (!quit) {
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
      send_message(connection, "WELCOME to Reindeer socket server\n", data.verbose);

      // start listen on multiple connection
      bool endconnection = false;

      while (!endconnection) {
          message = wait_message(connection, data.verbose);

              // analysis args coming from socket message
              switch(getCmd(message)) {
                  case Command::QUIT:
                      // QUIT message: Send a message to the connection and Quit
                      send_message(connection, "I'm leaving, see you next time !", data.verbose);
                      quit = true;
                      break;
                  case Command::STOP:
                      // stop client connection
                      send_message(connection, "See you soon !", data.verbose);
                      endconnection = true;
                      break;
                  case Command::HELP:
                      // HELP message
                      send_message(connection, help_command, data.verbose);
                      break;
                  case Command::INDEX:
                      // INDEX message: ask which index is in memory
                      message = "INDEX:";
                      message.append(reindeer_index.matrix_name);
                      send_message(connection, message, data.verbose);
                      break;
                  case Command::FILECMD:
                  // informations are by pair
                  { // to be able to initialize subcmds var, case are like goto
                      string query_file {};
                      string outFile {};
                      vector<string> subcmds = split_utils(message, ':');
                      if (not (subcmds.size() % 2)) {
                          if (data.verbose)
                              cerr << "in File mode, with paired arguments" << endl;
                          for (auto it = subcmds.begin(); it < subcmds.end(); it++) {
                              switch(getCmd(string(*it))) {
                                  case Command::FILECMD:
                                      // get fasta file
                                      it++;
                                      query_file = *it;
                                      cerr << "FILE = " << query_file << endl;
                                      break;
                                  case Command::THRESHOLD:
                                      // get threshold
                                      it++;
                                      reindeer_index.threshold = stoi(*it);
                                      cerr << "THRESHOLD = " << reindeer_index.threshold  << endl;
                                      break;
                                  case Command::OUTFILE:
                                      // get output file
                                      it++;
                                      outFile = *it;
                                      cerr << "OUTFILE = " << outFile << endl;
                                      break;
                                  case Command::FORMAT:
                                      // get output format
                                      it++;
                                      reindeer_index.output_format = *it;
                                      cerr << "FORMAT = " << reindeer_index.output_format << endl;
                                      break;
                                  default :
                                      message.insert(0, "UNKNOWN command: ");
                                      send_message(connection, message, data.verbose);
                              }
                          }
                          if (fs::exists(query_file)) {
                              if (!outFile.empty()) reindeer_index.output = outFile;
                              reindeer_index.querying(query_file, reindeer_index.threshold, reindeer_index.output_format);
                              // we have finished
                              send_message(connection, "DONE", data.verbose);
                          } else {
                              send_message(connection, "ERROR:The entry is not a file or not given (FILE:path.fa is required)", data.verbose);
                          }
                      } else {
                          send_message(connection, "ERROR: no file given to command FILE, or not paired command", data.verbose);
                      }  // if subcmds.size %2
                      break;
                  }
                  default :
                      message.insert(0, "UNKNOWN command: ");
                      send_message(connection, message, data.verbose);
              }
      } //while endconnection
      // Close the connections
      close(connection);
      } // while quit
    } // end of server mode
    close(sockfd);
    return 0;
}
