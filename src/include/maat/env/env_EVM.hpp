#ifndef MAAT_ENV_EVM_HPP
#define MAAT_ENV_EVM_HPP

#include "maat/env/library.hpp"
#include "maat/env/filesystem.hpp"
#include "maat/env/os.hpp"
#include "maat/env/env.hpp"
#include "maat/arch.hpp"
#include "maat/snapshot.hpp"
#include "maat/process.hpp"
#include "maat/serializer.hpp"
#include <iostream>

namespace maat{
namespace env{
/// Simulation of the Ethereum blockchain state
namespace EVM{


/** \addtogroup env
 * \{ */

/// EVM Stack
class Stack
{
private:
    std::vector<Value> _stack;
public:
    Stack() = default;
    virtual ~Stack() = default;
public:
    /// Return number of elements present in the stack
    int size() const;
    /// Get value at given position. Position 0 is the top of the stack
    const Value& get(int pos) const;
    /// Remove 'n' values from top of the stack
    void pop(int n=1);
    /// Get value at given position. Position 0 is the top of the stack
    void set(const Value& value, int pos);
    /// Push new value at the top of the stack
    void push(const Value& value);
public:
    friend std::ostream& operator<<(std::ostream&, const Stack&);
private:
    // Convert a position to corresponding index in internal vector
    int _pos_to_idx(int pos) const;
};

/// EVM Volatile Memory
class Memory
{
private:
    MemEngine _mem;
    addr_t _size; // Current memory size in the EVM sense
    addr_t _limit; // Limit of internally allocated memory
    addr_t _alloc_size;
    std::shared_ptr<VarContext> _varctx;
public:
    Memory(std::shared_ptr<VarContext> ctx);
    ~Memory() = default;
public:
    MemEngine& mem(); ///< Get internal memory engine
    addr_t size() const; ///< Get current memory size in bytes
public:
    /// Read 'nb_bytes' at address 'addr'
    Value read(const Value& addr, size_t nb_bytes);
    /// Write value to memory
    void write(const Value& addr, const Value& val);
public:
    /// Expand memory if needed to write 'nb_bytes' at 'addr'
    void expand_if_needed(const Value& addr, size_t nb_bytes);
};


class ValueHash {
  public:
    ::std::size_t operator ()(const Value& value) const
    {
        if (value.is_abstract())
            return value.as_expr()->hash();
        else
            return value.as_uint(); // Not abstract so should be safe to call
    }
};

class ValueEqual {
  public:
    bool operator ()(const Value& v1, const Value& v2) const
    {
        return v1.eq(v2);
    }
};

/// Contract permananent storage
class Storage
{
private:
    /// Storage state for concrete addresses
    std::unordered_map<Value, Value, ValueHash, ValueEqual> _storage;
    /// History of storage writes, including symbolic addresses
    std::vector<std::pair<Value, Value>> writes_history;
    std::shared_ptr<VarContext> _varctx;
    /// True if at least one address written to was symbolic
    bool _has_symbolic_addresses; 
public:
    Storage(std::shared_ptr<VarContext> ctx);
    ~Storage() = default;
public:
    /// Get storage word at 'addr'
    Value read(const Value& addr);
    /// Write storage word at 'addr'
    void write(const Value& addr, const Value& val, const Settings& settings);
};

/// Result of a call inside a contract
class TransactionResult
{
public:
    enum class Type : uint8_t 
    {
        RETURN,
        REVERT,
        STOP,
        NONE
    };
private:
    Type _type;
    std::vector<Value> _return_data;
public:
    TransactionResult(Type type, std::vector<Value> return_data);
public:
    const std::vector<Value>& return_data() const;
    size_t return_data_size() const;
    Type type() const;
};

/// Call into a smart contract
class Transaction
{
public:
    Value origin; ///< Original account that initiated contract execution
    Value sender; ///< Account/contract that initiated this transaction
    Number recipient; ///< Recipient of the transaction
    Value value; ///< Number of ether to transfer from sender to recipient
    std::vector<Value> data; ///< Additionnal transaction data
    Value gas_limit; ///< Maximum amount of gas that can be consumed by the transaction
    std::optional<TransactionResult> result; ///< Result of the transaction
public:
    Transaction();
    Transaction(
        Value origin,
        Value sender,
        Number recipient,
        Value value,
        std::vector<Value> data,
        Value gas_limit
    );
public:
    /// Total size of transaction data in bytes
    size_t data_size() const;
    /// Return 32 bytes from transaction data starting at byte 'offset'
    Value data_load_word(size_t offset) const;
    /// Return bytes [offset ... offset+len] from transaction data
    std::vector<Value> data_load_bytes(size_t offset, size_t len) const;
};

/// Deployed Smart-Contract
class Contract
{
public:
    Value address; ///< Address where the contract is deployed
    Stack stack; ///< Stack of the executing EVM
    Memory memory; ///< Volatile memory of the executing EVM
    Storage storage; ///< Persistent contract storage
    std::optional<Transaction> transaction; ///< Transaction being executed
protected:
    std::optional<TransactionResult> result_from_last_call; ///< Result of last call emitted from the current executing environment
public:
    unsigned int code_size; ///< Size of code currently executing
public:
    /** Constructor. Create new deployed contract */
    Contract(const MaatEngine& engine, Value address);
};

// Util functions
void _set_EVM_code(MaatEngine& engine, const std::vector<Value>& code);
void _set_EVM_code(MaatEngine& engine, uint8_t* code, size_t code_size);

typedef std::shared_ptr<Contract> contract_t;


/// Helper class for simulating the KECCAK hash function symbolically
class KeccakHelper
{
private:
    std::string _symbolic_hash_prefix;
    std::unordered_map<Value, Value, ValueHash, ValueEqual> known_hashes;
public:
    KeccakHelper();
public:
    /// Return the result of applying the KECCAK hash function to 'val'
    Value apply(VarContext& ctx, const Value& val, uint8_t* raw_bytes);
    /// Get the prefix used for symbolic hash results
    const std::string& symbolic_hash_prefix() const;
};

// Compute the keccak hash of bytes and return the result as a 'Value'
Value _do_keccak256(uint8_t* in, int size);

/// Specialisation of 'EnvEmulator' for the Ethereum blockchain state
class EthereumEmulator: public EnvEmulator
{
private:
    int _uid_cnt;
    std::unordered_map<int, contract_t> _contracts;
public:
    KeccakHelper keccak_helper;
public:
    EthereumEmulator();
    virtual ~EthereumEmulator() = default;
private:
    /// In-place initialization function used by constructor and deserializer
    void _init();
public:
    /// Add a running deployed contract instance and return its unique id
    int add_contract(contract_t contract);
    /// Get running contract by uid
    contract_t get_contract_by_uid(int uid) const;
public:
    virtual maat::serial::uid_t class_uid() const;
    virtual void dump(maat::serial::Serializer& s) const;
    virtual void load(maat::serial::Deserializer& d);
};

/** \brief Helper function that gets the environment linked to an engine and casts it to EthereumEmulator.
 * This function performs dynamic casting without further checks, use at your own risk */ 
std::shared_ptr<EthereumEmulator> get_ethereum(MaatEngine& engine);

/// Helper function that gets the running contract associated to an engine
contract_t get_contract_for_engine(MaatEngine& engine); 

/// Helper function that converts a hex string into bytes. Eg "0102" -> "\x01\x02"
std::vector<uint8_t> hex_string_to_bytes(const std::vector<char>& in);

/** \} */ // doxygen group env

} // namespace EVM
} // namespace env
} // namespace maat

#endif 