#include "maat/serializer.hpp"
#include "maat/maat.hpp"

namespace maat{
namespace serial{

Deserializer::Deserializer(std::istream& in): _stream(in){}

Deserializer::Stream& Deserializer::stream() {return _stream;}

void Deserializer::init()
{
    // Read index info
    int index_pos, index_cnt;
    uid_t root_obj_uid;
    stream() >> bits(index_pos) >> bits(index_cnt) >> bits(root_obj_uid);
    uid_to_object.reserve(index_cnt);
    data_pos_to_object.reserve(index_cnt);

    // Go through every index entry and create empty object
    stream().set_pos(index_pos);
    for (int i = 0; i < index_cnt; i++)
    {
        Serializer::IndexEntry entry{0,0,0,0};
        stream() >> bits(entry.obj_uid) >> bits(entry.class_uid)
                 >> bits(entry.data_pos) >> bits(entry.data_end_pos);
        // Save data position of root object
        if (entry.obj_uid == root_obj_uid)
        {
            root_obj_data_pos = entry.data_pos;
            root_obj_data_end_pos = entry.data_end_pos;
        }
        // Make sure that uid is unique
        if (uid_to_object.find(entry.obj_uid) != uid_to_object.end())
            throw serialize_exception("Deserializer::init(): got non unique uid!");
        // Create object
        Serializable* obj_ptr = _factory.new_object(entry.class_uid);
        // Associate uid to real object
        uid_to_object[entry.obj_uid] = obj_ptr;
        data_pos_to_object[entry.data_pos] = obj_ptr;
    }
}

Serializable* Deserializer::_deserialize(bool skip_root_obj)
{
    int index_pos = 0, index_cnt = 0;
    uid_t root_obj_uid = 0;
    // Initialise the objects
    init();
    // Set stream to beginning of data (skip index_pos, index_cnt)
    stream().set_pos(0);
    stream() >> bits(index_pos) >> bits(index_cnt) >> bits(root_obj_uid);
    // Load all object contents
    try
    {
        while (stream().current_pos() < index_pos)
        {
            if (skip_root_obj && stream().current_pos() == root_obj_data_pos)
            {
                stream().set_pos(root_obj_data_end_pos);
            }
            else
            {
                Serializable* obj = data_pos_to_object.at(stream().current_pos());
                obj->load(*this);
            }
        }
    }
    catch (const std::out_of_range&)
    {
        throw serialize_exception("Deserializer::deserialize(): data position in stream doesn't correspond to an object");
    }

    // If we skipped root object, delete it from the map, because we won't use
    // it to create and return a pointer.
    if (skip_root_obj)
        delete uid_to_object[root_obj_uid];

    // TODO(boyan): check that every raw pointer now has an owner 
    // (either unique_ptr or shared_ptr)

    // Return root object
    return uid_to_object.at(root_obj_uid);
}


void Deserializer::deserialize(Serializable& dest)
{
    _deserialize(true);
    stream().set_pos(root_obj_data_pos);
    dest.load(*this);
}

Deserializer& Deserializer::operator>>(std::string& str)
{
    size_t size = 0;
    stream() >> bits(size); // read size of string
    std::vector<char> contents(size);
    stream() >> buffer(contents.data(), size); // read content
    str.assign(contents.data(), size); // set string
    return *this;
}

Serializable* Deserializer::Factory::new_object(uid_t class_uid)
{
    switch (class_uid)
    {
        case ClassId::ARCH_NONE:
            return new ArchNone();
        case ClassId::ARCH_X64:
            return new X64::ArchX64();
        case ClassId::ARCH_X86:
            return new X86::ArchX86();
        case ClassId::CONSTRAINT:
            return new ConstraintObject();
        case ClassId::CPU:
            return new ir::CPU(0);
        case ClassId::ENV_EMULATOR:
            return new env::EnvEmulator();
        case ClassId::ENV_ETHEREUM_EMULATOR:
            return new env::EVM::EthereumEmulator();
        case ClassId::ENV_LINUX_EMULATOR:
            return new env::LinuxEmulator(Arch::Type::NONE);
        case ClassId::EVM_CONTRACT:
            return new env::EVM::Contract();
        case ClassId::EXPR_BINOP:
            return new ExprBinop();
        case ClassId::EXPR_CONCAT:
            return new ExprConcat();
        case ClassId::EXPR_CST:
            return new ExprCst();
        case ClassId::EXPR_EXTRACT:
            return new ExprExtract();
        case ClassId::EXPR_ITE:
            return new ExprITE();
        case ClassId::EXPR_UNOP:
            return new ExprUnop();
        case ClassId::EXPR_VAR:
            return new ExprVar();
        case ClassId::FS_DIRECTORY:
            return new env::Directory(nullptr);
        case ClassId::INTERVAL_TREE:
            return new IntervalTree();
        case ClassId::LIFTER:
            return new Lifter(CPUMode::NONE);
        case ClassId::MAAT_ENGINE:
            return new MaatEngine(Arch::Type::NONE);
        case ClassId::MEM_ABSTRACT_BUFFER:
            return new MemAbstractBuffer();
        case ClassId::MEM_CONCRETE_BUFFER:
            return new MemConcreteBuffer();
        case ClassId::MEM_ENGINE:
            return new MemEngine();
        case ClassId::MEM_SEGMENT:
            return new MemSegment(0, 0);
        case ClassId::MEM_STATUS_BITMAP:
            return new MemStatusBitmap();
        case ClassId::PATH_MANAGER:
            return new PathManager();
        case ClassId::PHYSICAL_FILE:
            return new env::PhysicalFile(nullptr, env::PhysicalFile::Type::INVALID);
        case ClassId::PROCESS_INFO:
            return new ProcessInfo();
        case ClassId::SNAPSHOT_MANAGER:
            return new SnapshotManager<maat::Snapshot>();
        case ClassId::SNAPSHOT_MANAGER_ENV:
            return new SnapshotManager<maat::env::Snapshot>();
        case ClassId::SYMBOL_MANAGER:
            return new SymbolManager();
        case ClassId::SYMBOLIC_MEM_ENGINE:
            return new SymbolicMemEngine(0, nullptr, Endian::LITTLE);
        case ClassId::VALUE:
            return new Value();
        case ClassId::VAR_CONTEXT:
            return new VarContext();
        default:
            throw serialize_exception(
                Fmt() << "Deserializer::Factory::new_object: unsupported class uid: " 
                      << (int)class_uid >> Fmt::to_str
            );
    }
}

std::shared_ptr<Serializable> Deserializer::Factory::new_shared_ptr(Serializable* raw_ptr)
{
    if (already_has_unique_ptr.find(raw_ptr) != already_has_unique_ptr.end())
        throw serialize_exception(
            "Deserializer: Trying to create shared_ptr from raw pointer that was already used to create a unique_ptr"
        );

    auto it = obj_to_shared_ptr.find(raw_ptr);
    if (it != obj_to_shared_ptr.end())
        return it->second;
    else
    {
        obj_to_shared_ptr[raw_ptr] = std::shared_ptr<Serializable>(raw_ptr);
        return obj_to_shared_ptr[raw_ptr];
    }
}

std::unique_ptr<Serializable> Deserializer::Factory::new_unique_ptr(Serializable* raw_ptr)
{
    if (already_has_unique_ptr.find(raw_ptr) != already_has_unique_ptr.end())
        throw serialize_exception(
            "Deserializer: Trying to create unique_ptr from raw pointer that was already used to create another unique_ptr"
        );

    auto it = obj_to_shared_ptr.find(raw_ptr);
    if (it != obj_to_shared_ptr.end())
        throw serialize_exception(
            "Deserializer: Trying to create unique_ptr from raw pointer that was already used to create a shared_ptr"
        );
    else
    {
        already_has_unique_ptr.insert(raw_ptr);
        return std::unique_ptr<Serializable>(raw_ptr);
    }
}

Deserializer::Stream::Stream(std::istream& _in): in(std::ref(_in)) {}

int Deserializer::Stream::current_pos() const
{
    return (int)in.get().tellg();
}

void Deserializer::Stream::set_pos(int pos)
{
    in.get().seekg((std::streampos)pos);
}

} // namespace serial
} // namespace maat