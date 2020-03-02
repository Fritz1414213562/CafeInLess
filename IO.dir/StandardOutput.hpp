#ifndef STANDARD_OUTPUT_HPP
#define STANDARD_OUTPUT_HPP
#include<iostream>
#include<type_traits>
#include<string>
#include<tuple>

// dump to standard output

namespace CafeInLess {

namespace IO {

class Standard_Output {

public:
    Standard_Output() = default;
    ~Standard_Output() = default;

    template<typename... Args>
    void operator()(Args... args) {
        output_Messages(args...);
    }

//	void operator[](const std::string& message, const std::size_t& hyphen_num) {
//		output_HyphenBlock(message, hyphen_num);
//	}

	void operator[](const std::size_t& hyphen_num) {
		output_HyphenBlock("", hyphen_num);
	}


private:
    

    template<typename message_type>
    void output_Message(const message_type& message) {
        std::cout << message;
    }


    template<typename... Types>
    void output_Messages(Types... m_messages) {
        constexpr std::size_t Types_size = sizeof...(Types);
        static_assert(Types_size > 0, "No arguments");
        std::tuple<Types...> message_tuple = std::make_tuple(m_messages...);
        iterate_Message(message_tuple);
    }


    template<std::size_t iterative_num = 0, typename Tuple_Type>
    void iterate_Message(const Tuple_Type& Tuple) {
        if constexpr (iterative_num < std::tuple_size<Tuple_Type>::value - 1) {
            output_Message(std::get<iterative_num>(Tuple));
			std::cout << " ";
            iterate_Message<iterative_num + 1>(Tuple);
        }
		else if constexpr (iterative_num == std::tuple_size<Tuple_Type>::value - 1) {
			output_Message(std::get<iterative_num>(Tuple));
			std::cout << std::endl;
		}
    }


    void output_HyphenBlock(const std::string& message, const std::size_t hyphen_num) {
        const std::size_t iterate_num = (message.size() < hyphen_num) * hyphen_num + (message.size() >= hyphen_num) * message.size() - message.size();
        std::cout << "------" << message;
        for (std::size_t idx = 0; idx < iterate_num; ++idx) {
            std::cout << "-";
        }
        std::cout << std::endl;
    }

    

};
}
}

#endif /* STANDARD_OUTPUT_HPP */
