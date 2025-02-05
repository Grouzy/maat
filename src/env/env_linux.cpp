#include "maat/env/env.hpp"

namespace maat{
namespace env{

LinuxEmulator::LinuxEmulator(Arch::Type arch): EnvEmulator(arch, OS::LINUX)
{
    _init(arch);
}

void LinuxEmulator::_init(Arch::Type arch)
{
    EnvEmulator::_init(arch, OS::LINUX);
    _arch = arch;
    // Load emulated libraries
    switch (arch)
    {
        case Arch::Type::X86:
            _libraries.push_back(env::emulated::linux_x86_libc());
            _syscall_func_map = env::emulated::linux_x86_syscall_map();
            break;
        case Arch::Type::X64:
            _libraries.push_back(env::emulated::linux_x64_libc());
            _syscall_func_map = env::emulated::linux_x64_syscall_map();
            break;
        case Arch::Type::NONE:
        default:
            break;
    }
}

void LinuxEmulator::add_running_process(const ProcessInfo& pinfo, const std::string& filepath)
{
    // Create actual file
    fs.add_real_file(filepath, pinfo.binary_path, true);

    // Set symbolic links to loaded binary in /proc/<pid>/exe
    // and /proc/self/exe
    std::stringstream ss;
    ss << "/proc/" << std::dec << pinfo.pid << "/exe";
    fs.create_symlink(ss.str(), pinfo.binary_path, true);
    fs.create_symlink("/proc/self/exe", pinfo.binary_path, true);

    // Create stdin,stdout,stderr, for this process
    std::string _stdin = fs.get_stdin_for_pid(pinfo.pid);
    std::string _stdout = fs.get_stdout_for_pid(pinfo.pid);
    std::string _stderr = fs.get_stderr_for_pid(pinfo.pid);

    fs.create_file(_stdin);
    fs.create_file(_stdout);
    fs.create_file(_stderr);
    fs._new_fa(_stdin, 0);
    fs._new_fa(_stdout, 1);
    fs._new_fa(_stderr, 2);
    fs.get_file_by_handle(1)->flush_stream = std::ref(std::cout);
    fs.get_file_by_handle(2)->flush_stream = std::ref(std::cerr);
}

uid_t LinuxEmulator::class_uid() const
{
    return serial::ClassId::ENV_LINUX_EMULATOR;
}

void LinuxEmulator::dump(serial::Serializer& s) const
{
    s << bits(_arch) << fs;
}

void LinuxEmulator::load(serial::Deserializer& d)
{
    Arch::Type arch;
    d >> bits(arch);
    _init(arch);
    // Then after initialisation, load filesystem
    d >> fs;
}

} // namespace env
} // namespace maat
