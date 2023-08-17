#include <sys/socket.h> // For socket functions
#include <netinet/in.h> // For sockaddr_in
#include <getopt.h>
#include <unistd.h>
#include "reindeer.hpp"
#include "utils.hpp"

using namespace std;
namespace fs = filesystem;

typedef struct data_for_socket{
    int port = 0;
    string index_directory = "";
} DataSoc;

void help() {
    cout << endl;
    cout <<
    "* Mandatory parameters *\n"
    "-p, --port <number>                :       Port\n"
    "-l <directory>                     :       Index directory\n"
    "\n"
    "* General parameters *\n"
    "-h, --help                         :       Print help (what you are currently seeing)\n"
    << endl;
    exit(0);
}

void process_args(int argc, char** argv, DataSoc& data) {

    const char* const short_opts = "hp:l:";
    const option long_opts[] = {
        { "port", required_argument, nullptr, 'p' },
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
            case 'l':
                data.index_directory = optarg;
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
        cout << "Failed to create socket. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }

    // avoid err 98: socket already in use, due to timeout when server close
    int yes = 1;
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (void*)&yes, sizeof(yes)) < 0) {
        cout << "setsockopt() failed. Error: errno: " << errno << endl;
    }
    // Listen to port data.port on any address
    sockaddr_in sockaddr;
    sockaddr.sin_family = AF_INET;
    sockaddr.sin_addr.s_addr = INADDR_ANY;
    sockaddr.sin_port = htons(data.port); // htons is necessary to convert a number to
    // network byte order
    if (bind(sockfd, (struct sockaddr*)&sockaddr, sizeof(sockaddr)) < 0) {
        cout << "Failed to bind to port " << data.port << ". errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }

    // Start listening. Hold at most 10 connections in the queue
    if (listen(sockfd, 1) < 0) {
        cout << "Failed to listen on socket. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }

    cout << " * created socket on port " << data.port << endl;
    cout << " * waiting for client to connect" << endl;
    // Grab a connection from the queue
    auto addrlen = sizeof(sockaddr);
    int connection = accept(sockfd, (struct sockaddr*)&sockaddr, (socklen_t*)&addrlen);
    if (connection < 0) {
        cout << "Failed to grab connection. errno: " << errno << endl;
        exit(EXIT_FAILURE);
    }
    cout << " * connected to one client" << endl;
    // give a welcome message
    message = "WELCOME to Reindeer socket server\n";
    write(connection, message.c_str(), message.size());

    string output = "reindeer_index_files", query_file = "", format = "raw";
    // create reindeer index object
    Reindeer_Index<uint16_t> reindeer_index(data.index_directory, output, 40, query_file, 1, true, format);
    //load index
    reindeer_index.load_index();

    // give which index is used
    message = "INDEX:";
    message.append(reindeer_index.matrix_name);
    message.append("\n");
    write(connection, message.c_str(), message.size());

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
                    message = "See you soon !";
                    cout << message << endl;
                    write(connection, message.c_str(), message.size());
                    quit = true;
                    break;
                // HELP message
                case 'H':
                    message = " * HELP = help message\n";
                    message.append(" * QUIT = quit message\n");
                    message.append(" * INDEX = ask index in use\n");
                    message.append(" * FILE:myfile.fasta[:THRESHOLD:value][:OUTFILE:myoutfile]\n");
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
                    vector<string> subcmds = split_utils(strBuffer, ':');
                    if (not (subcmds.size() % 2)) {
                        for (auto it = subcmds.begin(); it < subcmds.end(); it++) {
                            char option = toupper(string(*it).at(0));
                            switch(option) {
                                // get fasta file: FILE
                                case 'F':
                                    it++;
                                    reindeer_index.query = *it;
                                    cout << "FILE = " << reindeer_index.query << endl;
                                    break;
                                // get threshold: THRESHOLD
                                case 'T':
                                    it++;
                                    reindeer_index.threshold = stoi(*it);
                                    cout << "THRESHOLD = " << reindeer_index.threshold  << endl;
                                    break;
                                // get output file: OUTFILE
                                case 'O':
                                    it++;
                                    outFile = *it;
                                    cout << "OUTFILE = " << outFile << endl;
                                    break;
                                default :
                                    message = "UNKNOWN command: ";
                                    message.append(1, option);
                                    message.append("\n");
                                    cout << message << endl;
                                    write(connection, message.c_str(), message.size());
                            }
                        }
                        if (fs::exists(reindeer_index.query)) {
                            if (!outFile.empty()) reindeer_index.output = outFile;
                            reindeer_index.querying();
                            // we have finished
                            message = "DONE\n";
                            write(connection, message.c_str(), message.length());
                        } else {
                            message = "ERROR:The entry is not a file or not given (FILE:path.fa is required)\n";
                            cout << message << endl;
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
    close(sockfd);
    return 0;
}
