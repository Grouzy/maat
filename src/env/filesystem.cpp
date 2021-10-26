#include "env/filesystem.hpp"
#include <fstream>

namespace maat
{
namespace env
{
   
// ================ PhysicalFile =====================
PhysicalFile::PhysicalFile(PhysicalFile::Type t): type(t)
{
    data = std::make_shared<MemSegment>(0x0, 0xfff); // Default size: 0x1000
    deleted = false;
    _size = 0;
    istream_read_offset = 0;
}

unsigned int PhysicalFile::size()
{
    return _size;
}

unsigned int PhysicalFile::copy_real_file(const std::string& filename)
{
    // Read file content
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<char> buffer(size);
    if (file.read(buffer.data(), size))
    {
        addr_t offset = 0;
        return write_buffer((uint8_t*)buffer.data(), offset, size);
    }
    else
    {
        throw env_exception(
            Fmt() << "Error reading contents of '" << filename << "'"
            >> Fmt::to_str
        );
    }
}

unsigned int PhysicalFile::write_buffer(const std::vector<Expr>& buffer, addr_t& write_ptr)
{
    int n = 0;
    VarContext dummy_ctx;

    _adjust_write_offset(write_ptr);

    if (deleted)
    {
        throw env_exception("Trying to write to deleted file");
    }

    if (is_symlink())
    {
        throw env_exception("Can not write to symbolic link file");
    }
    
    for (const auto& e : buffer)
    {
        if (e->size/8 + write_ptr - 1 > data->end)
        {
            // Extend file to write more
            data->extend_after(data->size()); // Double size
        }
        /* TODO
        if( _snapshot_manager != nullptr ){
            // If snapshot manager provided, record write
            _snapshot_manager->record_fs_file_write(this, state.write_ptr, e->size/8);
        } */
        data->write(write_ptr, e, dummy_ctx);
        write_ptr += e->size/8;
        n += e->size/8;
    }

    // Update size
    if (write_ptr > _size)
    {
        _size = write_ptr;
    }

    return n;
}

unsigned int PhysicalFile::write_buffer(uint8_t* buffer, addr_t& offset, int nb_bytes)
{
    if (deleted)
    {
        throw env_exception("Trying to write to deleted file");
    }

    if (is_symlink())
    {
        throw env_exception("Can not write to symbolic link file");
    }
    
    if (nb_bytes == 0)
        return 0;

    _adjust_write_offset(offset);

    while (offset + nb_bytes -1 > data->end)
    {
        // Extend file to write more
        data->extend_after(nb_bytes > data->size()? nb_bytes : data->size());
    }

    /* TODO
    if( _snapshot_manager != nullptr ){
        // If snapshot manager provided, record write
        _snapshot_manager->record_fs_file_write(this, offset, nb_bytes);
    } */

    data->write(offset, buffer, nb_bytes);

    // Update size
    offset += nb_bytes;
    if (offset > _size)
    {
        _size = offset;
    }

    return nb_bytes;
}

unsigned int PhysicalFile::read_buffer(
    std::vector<Expr>& buffer,
    addr_t& read_ptr,
    unsigned int nb_elems,
    unsigned int elem_size
)
{
    unsigned int cnt = 0;
    if (deleted)
    {
        throw env_exception("Trying to read from deleted file");
    }

    if (is_symlink())
    {
        throw env_exception("Can not read from symbolic link file");
    }

    _adjust_read_offset(read_ptr);

    // If nothing to read, just return
    if (read_ptr >= _size)
    {
        return 0;
    }
    // Else do read
    for (int i = 0; i < nb_elems; i++)
    {
        if (read_ptr + elem_size > _size)
        {
            // Reached end of file buffer
            buffer.push_back(data->read(read_ptr, _size-read_ptr));
            cnt += (_size-read_ptr);
            read_ptr = _size;
            break;
        }
        else
        {
            // Normal read
            buffer.push_back(data->read(read_ptr, elem_size));
            read_ptr += elem_size;
            cnt += elem_size;
        }
    }
    // Update the read offset (from streams)
    istream_read_offset = read_ptr;

    return cnt;
}

void PhysicalFile::_adjust_write_offset(addr_t& offset)
{
    if (type == PhysicalFile::Type::IOSTREAM)
        offset = _size;
}

void PhysicalFile::_adjust_read_offset(addr_t& offset)
{
    if (type == PhysicalFile::Type::IOSTREAM)
        offset = istream_read_offset;
}


bool PhysicalFile::is_symlink()
{
    return type == PhysicalFile::Type::SYMLINK;
}

const std::string& PhysicalFile::symlink()
{
    if (not is_symlink())
        throw env_exception("PhysicalFile::symlink(): File is not a symlink!");
    return _symlink;
}

void PhysicalFile::_set_symlink(const std::string& target)
{
    type = PhysicalFile::Type::SYMLINK;
    _symlink = target;
}

bool PhysicalFile::is_deleted()
{
    return deleted;
}

void PhysicalFile::set_deleted(bool val)
{
    deleted = val;
}

FileAccessor::FileAccessor(physical_file_t file, filehandle_t handle):
flags(0), physical_file(file), _handle(handle), deleted(false), _alloc_addr(0)
{
    // Init access state
    state.read_ptr = 0;
    state.write_ptr = 0;
}

unsigned int FileAccessor::write_buffer(const std::vector<Expr>& buffer)
{
    return physical_file->write_buffer(buffer, state.write_ptr);
}

unsigned int FileAccessor::write_buffer(uint8_t* buffer, int len)
{
    return physical_file->write_buffer(buffer, state.write_ptr, len);
}

unsigned int FileAccessor::read_buffer(
    std::vector<Expr>& buffer,
    unsigned int nb_elems,
    unsigned int elem_size
){
    return physical_file->read_buffer(buffer, state.read_ptr, nb_elems, elem_size);
}

filehandle_t FileAccessor::handle() const
{
    return _handle;
}

// ====================== Directory ======================
Directory::Directory(): deleted(false){};

bool Directory::create_file(fspath_t path, bool create_path)
{
    if (path.empty())
    {
        std::cout << "DEBUG EMPTY PATH " << std::endl;
        return false;
    }
    else if (path.size() == 1)
    {
        
        std::cout << "DEBUG path[-1] " << path.back() << std::endl;
        // Check if file name already taken
        if (_contains_name(path.back()))
        {
            std::cout << "DEBUG already exists " << std::endl;
            return false; // File already exists ! 
        }

        // Create new file
        physical_file_t file = std::make_shared<PhysicalFile>();
        files[path.back()] = file;
        return true;
    }
    else
    {
        std::cout << "DEBUG subdir " << path.front() << std::endl;
        // Create file in subdir
        for (auto& dir : subdirs)
        {
            if (dir.first == path.front())
            {
                path.erase(path.begin());
                return dir.second->create_file(path, create_path);
            }
        } 
        // No matching subdir, check if create_path is set
        if (create_path)
        {
            // Create dir
            std::string new_dir = path.front();
            if (!create_dir({new_dir}))
            {
                return false;
            }
            path.erase(path.begin());
            return get_dir({new_dir})->create_file(path, create_path);
        }
        // Fail :'(
        return false;
    }
}

bool Directory::create_dir(fspath_t path)
{
    if (path.empty())
    {
        return false;
    }
    else if (path.size() == 1)
    {
        // Check if name already taken
        if (_contains_name(path.back()))
            return false;

        // Create new dir
        directory_t dir = std::make_shared<Directory>();
        subdirs[path.back()] = dir;
        return true;
    }
    else
    {
        // Create dir in subdir
        for (auto& subdir : subdirs)
        {
            if (subdir.first == path.front())
            {
                path.erase(path.begin());
                return subdir.second->create_file(path);
            }
        }
        // No matching subdir :'(
        return false;
    }
}

physical_file_t Directory::get_file(fspath_t path)
{
    if (path.empty())
    {
        throw env_exception("File doesn't exist");
    }
    else if (path.size() == 1)
    {
        // Check if file exists
        for (auto& file : files)
        {
            if (file.first == path.back())
                return file.second;
        }
        throw env_exception(
            Fmt() << "File '" << path.back() << "' doesn't exist"
            >> Fmt::to_str
        );
    }
    else
    {
        // Get file in subdir
        std::string dir_name = path.front();
        path.erase(path.begin());
        return get_dir({dir_name})->get_file(path);
    }
}

directory_t Directory::get_dir(fspath_t path)
{
    if (path.empty())
    {
        throw env_exception("Didn't find directory");
    }
    else if (path.size() == 1)
    {
        // Check if dir exists
        for (auto& subdir : subdirs)
        {
            if (subdir.first == path.back())
                return subdir.second;
        }
        throw env_exception(
            Fmt() << "Directory '" << path.back() << "' doesn't exist"
            >> Fmt::to_str
        );
    }
    else
    {
        // Get dir in subdir
        std::string dir_name = path.front();
        path.erase(path.begin());
        return get_dir({dir_name})->get_dir(path);
    }
}

bool Directory::delete_file(fspath_t path, bool weak)
{
    if (path.empty())
    {
        return false;
    }
    else if (path.size() == 1)
    {
        // Check if file
        auto it = files.find(path.back());
        if (it != files.end())
        {
            if (weak)
                it->second->set_deleted(true);
            else
                files.erase(it);
            return true;
        }
        return false; // File doesn't exist
    }
    else
    {
        // Delete file in subdir
        std::string dir_name = path.front();
        path.erase(path.begin());
        return get_dir({dir_name})->delete_file(path);
    }
}

bool Directory::delete_dir(fspath_t path, bool weak)
{
    if (path.empty())
    {
        return false;
    }
    else if (path.size() == 1)
    {
        // Check if dir
        auto it = subdirs.find(path.back());
        if (it != subdirs.end())
        {
            if (weak)
                it->second->deleted = true;
            else
                subdirs.erase(it);
            return true;
        }
        return false; // File doesn't exist
    }
    else
    {
        // Delete dir in subdir
        std::string dir_name = path.front();
        path.erase(path.begin());
        return get_dir({dir_name})->delete_dir(path);
    }
}


// TODO not really useful since delete is now a erase() or clear()...
void Directory::delete_self(bool recursive, bool weak)
{
    fspath_t tmp_path;
    if (not weak)
    {
        files.clear();
        subdirs.clear();
        return;
    }

    // Weak delete
    // Delete files
    for (auto& file : files)
    {
        if (not file.second->is_deleted())
        {
            file.second->set_deleted(true);
            /* TODO
            if( sm ){
                tmp_path = self_path; 
                tmp_path.push_back(file->name);
                sm->record_fs_action(FileSystemAction(ENV_FS_DELETE_FILE, fs->get_filename_from_path(tmp_path)));
            } */
        }
    }
    // Delete subdirs
    /* TODO 
    for (auto& dir : subdirs)
    {
        
        if (not dir.second->deleted){
            tmp_path = self_path; 
            tmp_path.push_back(dir->name);
            dir->_delete_self_recursive(tmp_path, fs, sm);
        }
    } */

    // Deleted self
    deleted = true;

    /* TODO 
    if( sm ){
        sm->record_fs_action(FileSystemAction(ENV_FS_DELETE_DIR, fs->get_filename_from_path(self_path)));
    } */
}


node_status_t Directory::get_node_status(fspath_t path)
{
    node_status_t status = 0;
    if (path.empty())
    {
        return node::none;
    }
    
    try // Try file
    {
        physical_file_t file = get_file(path);
        status |= node::exists;
        status |= node::is_file;
        if (file->is_symlink())
            status |= node::is_symlink;
        return status;
    }
    catch(const env_exception& e)
    {
        try // Try dir
        {
            directory_t dir = get_dir(path);
            status |= node::exists;
            status |= node::is_dir;
            return status;
        }
        catch(const env_exception& e)
        {
            return node::none;
        }
        
    }
    return node::none;
}

bool Directory::_contains_name(const std::string& name)
{
    return  files.find(name) != files.end()
            or subdirs.find(name) != subdirs.end();
}

void Directory::print(std::ostream& os, const std::string& indent_string) const
{
    std::string tmp_indent;
    std::string next_indent;
    // Start with subdirs  
    tmp_indent = indent_string + " \u2502  ";
    next_indent = indent_string + " \u2502  ";
    for (
        auto subdir = subdirs.begin();
        subdir != subdirs.end() and subdir != std::prev(subdirs.end());
        subdir++
    )
    {
        // Print dir
        os << indent_string << " \u251C" << "\u2500\u2500";
        os << " " << subdir->first << "/\n";
        subdir->second->print(os, next_indent);
    }
    // Print last subdir
    if (not subdirs.empty())
    {
        auto subdir = subdirs.rbegin();
        if (files.empty())
        {
            tmp_indent = indent_string + "    ";
            os << indent_string << " \u2514" << "\u2500\u2500";
        }
        else
        {
            tmp_indent = indent_string + " \u2502  ";
            os << indent_string << " \u251C" << "\u2500\u2500";
        }
        os << " " << subdir->first << "/\n";
        subdir->second->print(os, tmp_indent);
    }

    // Then print files in this dir
    for (
        auto file = files.begin();
        file != files.end() and file != std::prev(files.end());
        file++
    )
    {
        os << indent_string << " \u251C" << "\u2500\u2500";
        os << " " << file->first << std::endl;
    }
    // Print last file
    if (not files.empty())
    {
        os << indent_string << " \u2514" << "\u2500\u2500 " << files.rbegin()->first << std::endl;
    }
}


// =============== FileSystem =================
FileSystem::FileSystem(OS system): _handle_cnt(0), orphan_file_wildcard('#')
{
    switch (system)
    {
        case OS::LINUX:
        case OS::NONE:
            path_separator = "/";
            reserved_handles = {0,1,2}; // stdin, stdout, stderr
            break;
        case OS::WINDOWS:
            path_separator = "\\";
            break;
        default:
            throw runtime_exception("FileSystem constructor: unsupported OS!");
    }
}

bool FileSystem::create_file(const std::string& path, bool create_path)
{
    Directory& dir = (path[0] == orphan_file_wildcard) ? orphan_files : root;

    if (not dir.create_file(fspath_from_path(path), create_path))
    {
        return false;
    }
    /* TODO
    if( _snapshot_manager != nullptr ){
        // Record fot snapshots
        _snapshot_manager->record_fs_action(FileSystemAction(ENV_FS_CREATE_FILE, get_filename_from_path(path))); 
    } */
    return true;
}

physical_file_t FileSystem::get_file(const std::string& path, bool follow_symlink)
{
    Directory& dir = (path[0] == orphan_file_wildcard) ? orphan_files : root;
    physical_file_t res = dir.get_file(fspath_from_path(path));
    while (follow_symlink and res != nullptr and res->is_symlink())
    {
        res = get_file(res->symlink());
    }
    return res;
}
physical_file_t FileSystem::get_file_by_handle(filehandle_t handle)
{
    FileAccessor& fa = get_fa_by_handle(handle);
    return fa.physical_file;
}

bool FileSystem::file_exists(const std::string& path)
{
    node_status_t status = get_node_status(path);
    return node::check_is_file(status);
}

bool FileSystem::delete_file(const std::string& path, bool weak)
{
    Directory& dir = (path[0] == orphan_file_wildcard) ? orphan_files : root;

    if (not dir.delete_file(fspath_from_path(path), weak))
    {
        return false;
    }

    /* TODO 
    if( _snapshot_manager != nullptr ){
        // Record fot snapshots
        _snapshot_manager->record_fs_action(FileSystemAction(ENV_FS_DELETE_FILE, get_filename_from_path(path)));
    } */
    return true;
}

bool FileSystem::create_symlink(
    const std::string& link,
    const std::string& pointed_file,
    bool create_path
)
{
    create_file(link, create_path); // create_file() handles snapshoting
    physical_file_t file = get_file(link);
    file->_set_symlink(pointed_file);
    return true;
}

bool FileSystem::create_dir(const std::string& path)
{
    if (not root.create_dir(fspath_from_path(path)))
    {
        return false;
    }
    /* TODO
    if( _snapshot_manager != nullptr ){
        // Record fot snapshots
        _snapshot_manager->record_fs_action(FileSystemAction(ENV_FS_CREATE_DIR???, get_filename_from_path(path))); 
    } */
    return true;
}

directory_t FileSystem::get_dir(const std::string& path)
{
    directory_t res = root.get_dir(fspath_from_path(path));
    return res;
}

bool FileSystem::delete_dir(const std::string& path, bool weak)
{
    if (not root.delete_dir(fspath_from_path(path), weak))
    {
        return false;
    }

    /* TODO 
    if( _snapshot_manager != nullptr ){
        // Record fot snapshots
        _snapshot_manager->record_fs_action(FileSystemAction(ENV_FS_DELETE_DIR???, get_filename_from_path(path)));
    } */
    return true;
}


filehandle_t FileSystem::new_fa(const std::string& path)
{
    filehandle_t handle = get_free_handle();
    _new_fa(path, handle);
    return handle;
}

void FileSystem::_new_fa(const std::string& path, filehandle_t handle)
{
    physical_file_t file = get_file(path);
    fa_list.emplace_back(FileAccessor(file, handle));
}

FileAccessor& FileSystem::get_fa_by_handle(filehandle_t handle)
{
    for (auto& fa : fa_list)
        if (fa.handle() == handle)
            return fa;
    throw env_exception(
        Fmt() << "No file accessor with handle: " << (int)handle
        >> Fmt::to_str
    );
}

void FileSystem::delete_fa(filehandle_t handle, bool weak)
{
    if (weak)
    {
        get_fa_by_handle(handle).deleted = true;
        // TODO Snapshot ???
    }
    else
    {
        fa_list.remove_if(
            [handle](const FileAccessor& fa){return fa.handle() == handle;}
        );
    }
}


std::string FileSystem::path_from_fspath(const fspath_t& path)
{
    std::string res = "";
    for (const auto& s : path)
    {
        res += path_separator;
        res += s;
    }
    return res;
}

fspath_t FileSystem::fspath_from_path(const std::string& path)
{
    // Default from root node
    return fspath_from_path_relative_to(path, {});
}

std::string FileSystem::path_from_relative_path(std::string rel_path, std::string path_base)
{
    return path_from_fspath(
        fspath_from_path_relative_to(
            rel_path,
            fspath_from_path(path_base)
        )
    );
}

fspath_t FileSystem::fspath_from_path_relative_to(std::string rel_path, fspath_t path_base)
{
    std::vector<std::string> res;
    std::string s;
    bool empty;
    int i = 0, pos;
    
    // Orphan file
    if (rel_path[0] == orphan_file_wildcard)
    {
        rel_path.erase(0,1);
        return {rel_path};
    }

    // Parse filename for delimiter
    while ((pos = rel_path.find(path_separator)) != std::string::npos)
    {
        s = rel_path.substr(0, pos);

        // Check if string is empty 
        empty =  s.find_first_not_of(" \t\n\v\f\r") == std::string::npos;
        
        if( i == 0 && !empty)
            // Relative path
            res = path_base;

        // Add dir to path
        if( s == "." || empty)
        {
            // Pass
        }
        else if( s == ".." )
        {
            if( !res.empty())
                res.pop_back();
        }
        else
        {
            // Normal dir/file
            res.push_back(s);
        }

        rel_path.erase(0, pos + path_separator.length());
        i++;
    }

    if (rel_path.empty())
    {
        // No final filename after last delimiter
        throw env_exception(
            Fmt() << "FileSystem::fspath_from_path(): invalid filename '"
            << rel_path << "'" >> Fmt::to_str
        );
    }
    else
    {
        res.push_back(rel_path);
    }

    return res;
}

std::string FileSystem::pointed_path_from_symlink(std::string symlink_file)
{
    physical_file_t file = get_file(symlink_file);
    std::string res = symlink_file;
    while (file->is_symlink())
    {
        res = file->symlink();
        file = get_file(file->symlink());
    }
    return res;
}

bool FileSystem::is_relative_path(const std::string& path)
{
    return (
        path.substr(0, path_separator.size()) != path_separator
    );
}

node_status_t FileSystem::get_node_status(const std::string& path)
{
    Directory& dir = (path[0] == orphan_file_wildcard) ? orphan_files : root;
    return dir.get_node_status(fspath_from_path(path));
}

filehandle_t FileSystem::get_free_handle()
{
    while (std::find(
            reserved_handles.begin(),
            reserved_handles.end(),
            _handle_cnt) != reserved_handles.end())
    {
        _handle_cnt++;
    }
    return _handle_cnt++;
}

std::string FileSystem::get_stdin_for_pid(int pid)
{
    std::stringstream ss; 
    ss << orphan_file_wildcard << "stdin_pid_" << std::dec << pid;
    return ss.str();
}

std::string FileSystem::get_stdout_for_pid(int pid)
{
    std::stringstream ss; 
    ss << orphan_file_wildcard << "stdout_pid_" << std::dec << pid;
    return ss.str();
}

std::string FileSystem::get_stderr_for_pid(int pid)
{
    std::stringstream ss; 
    ss << orphan_file_wildcard << "stderr_pid_" << std::dec << pid;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const FileSystem& fs)
{
    os << "File system tree: \n";
    fs.root.print(os, "  ");

    os << "\nOther files: \n";
    fs.orphan_files.print(os, "  ");
    return os;
}


namespace node
{
 
bool check_is_file(node_status_t s)
{
    return s & node::is_file;
}

bool check_is_dir(node_status_t s)
{
    return s & node::is_dir;
}

bool check_is_symlink(node_status_t s)
{
    return s & node::is_symlink;
}

}


} // namespace env
} // namespace maat
