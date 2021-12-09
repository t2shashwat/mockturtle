/*!
  \file aig_algebraic_rewriting.hpp
  \brief AIG algebraric rewriting

  EPFL CS-472 2021 Final Project Option 1
*/

#pragma once

#include "../networks/aig.hpp"
#include "../views/depth_view.hpp"
#include "../views/topo_view.hpp"
#include "../views/fanout_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class aig_algebraic_rewriting_impl
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  typedef struct Node_sig {
    node n;
    signal s;
    bool cp;
  } Node_sig;
public:
  aig_algebraic_rewriting_impl( Ntk& ntk )
    : ntk( ntk )
  {
    static_assert( has_level_v<Ntk>, "Ntk does not implement depth interface." );
  }

  void run()
  {
    ntk.foreach_node( [&]( node n ){
      std::cout << " ----- "<< n << " -----" << std::endl;
      std::cout << "LEVE : " << ntk.level(n) << std::endl;
      ntk.foreach_fanin(n, [this]( auto const& f) { std::cout << ntk.get_node(f) << " signal :"<< ntk.is_complemented(f)<<std::endl;
                                                                          });
    });
    ntk.foreach_po([&](signal const &f){ std::cout<< "po signal : "<< ntk.get_node(f) << "  " << ntk.is_complemented(f)<<std::endl;});
    std::cout << " ----- ===================== " << std::endl;
    bool cont{true}; /* continue trying */
    while ( cont )
    {
      cont = false; /* break the loop if no updates can be made */
      ntk.foreach_gate( [&]( node n ){
      if ( try_algebraic_rules( n ) )
        {
          ntk.update_levels();
          cont = true;
        }
      });
    }
    ntk.foreach_node( [&]( node n ){
      std::cout << " ----- "<< n << " -----" << std::endl;
      ntk.foreach_fanin(n, [this]( auto const& f) { std::cout << ntk.get_node(f) << " signal :"<< ntk.is_complemented(f)<< std::endl;
                                                                          });
    });
    ntk.foreach_po([&](signal const &f){ std::cout<< "po signal after update: "<< ntk.get_node(f) << "  " <<ntk.is_complemented(f)<<std::endl;});
    

  }

private:
  /* Try various algebraic rules on node n. Return true if the network is updated. */
  bool try_algebraic_rules( node n )
  {
    if ( try_associativity( n ) )
      return true;
    if ( try_distributivity( n ) )
      return true;
    /* TODO: add more rules here... */

    return false;
  }
  std::array<Node_sig, 2> ordered_children( node const& n ) const
  {
    std::array<Node_sig, 2> children;
    ntk.foreach_fanin(n, [this, &children]( auto const& f, uint32_t i ) { children[i].n = ntk.get_node(f);
                                                                          children[i].s = f;
                                                                          children[i].cp = ntk.is_on_critical_path(children[i].n);});
    return children;
  }
    
  /* Try the associativity rule on node n. Return true if the network is updated. */
  bool try_associativity( node n )
  {
    std::cout << " NODE : "<< n << std::endl; 

    auto children = ordered_children(n);
    
    bool one_pi = (ntk.is_pi(children[0].n)) ^ (ntk.is_pi(children[1].n));//nothing to do if both children are pi or both nodes
    if(!one_pi){
      std::cout << "##### NO UPDATE ***BOTH PI OR NODE " << std::endl;
      return false;
    } 
    if(!ntk.is_pi(children[0].n) && ntk.is_complemented(children[0].s)){
      return false;
    }
    if(!ntk.is_pi(children[1].n) && ntk.is_complemented(children[1].s)){
      return false;
    }
 
    Node_sig c_on_critical_path, c_not_on_critical_path;
    if(children[0].cp == true){
      c_on_critical_path = children[0];
      c_not_on_critical_path = children[1];
    }
    else if(children[1].cp == true){
      c_on_critical_path = children[1];
      c_not_on_critical_path = children[0];
    }
    else {
      return false;
    }
   auto grand_children = ordered_children(c_on_critical_path.n);
   Node_sig gc_on_critical_path, gc_not_on_critical_path;
    if(grand_children[0].cp == true){
      gc_on_critical_path = grand_children[0];
      gc_not_on_critical_path = grand_children[1];
    }
    else if(grand_children[1].cp == true){
      gc_on_critical_path = grand_children[1];
      gc_not_on_critical_path = grand_children[0];
    }

    bool both_pi = (ntk.is_pi(grand_children[0].n)) & (ntk.is_pi(grand_children[1].n));// nothing to do if both children are pi
    
    if(both_pi){
      std::cout << "##### NO UPDATE ******** BOTH PI " << std::endl;
      return false;
    }
    if(!ntk.is_pi(grand_children[0].n) && ntk.is_complemented(grand_children[0].s)){
      return false;
    }
    if(!ntk.is_pi(grand_children[1].n) && ntk.is_complemented(grand_children[1].s)){
      return false;
    }
    signal short_ckt = gc_on_critical_path.s;

    auto hanging = ntk.create_and(gc_not_on_critical_path.s, c_not_on_critical_path.s);
    auto new_out = ntk.create_and(short_ckt, hanging);
    ntk.substitute_node(c_on_critical_path.n, hanging);   
    ntk.substitute_node(n, new_out);

    return true;
  }
  std::tuple<Node_sig, bool> check_common_grand_child(std::array<Node_sig, 2> grand_children0, std::array<Node_sig, 2> grand_children1)
  {
    bool polarity_check = true;
    for(int gr0 = 0; gr0 < 2; gr0++){
      for(int gr1 = 0; gr1 < 2; gr1++){
        if(grand_children0[gr0].n == grand_children1[gr1].n){
          polarity_check = ntk.is_complemented(grand_children0[gr0].s) ^ ntk.is_complemented(grand_children1[gr1].s);
          if(!polarity_check)
            return std::make_tuple(grand_children0[gr0], true);
        }
      }
    }
    return std::make_tuple(grand_children0[0], false);
  }
  std::array<Node_sig, 2> get_other_two_nodes(std::array<Node_sig, 2> grand_children0, std::array<Node_sig, 2> grand_children1, node common_gc)
  {
    std::array<Node_sig, 2> non_common_gc;
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        if(i == 0){
          if(grand_children0[j].n!=common_gc){
            non_common_gc[0] = grand_children0[j]; 
          }
        }
        if(i==1){
          if(grand_children1[j].n!=common_gc){
            non_common_gc[1] = grand_children1[j]; 
          }
        }
        
      }
    }
    return non_common_gc;
    
  }
  bool is_aig_or(node n)
  {
    auto children = ordered_children(n);
    node out_node;
    signal out_sig;
    bool out_sig_complemented = false;
    std::cout << " OUT SIG  090909009    " << ntk.depth()  << "  " << ntk.level(n)<<std::endl;
    if(ntk.depth() == ntk.level(n)){
       ntk.foreach_po([this, &out_sig](signal const& ns){ out_sig = ns;});
       if(ntk.is_complemented(out_sig) == true){
         out_sig_complemented = true;
         std::cout << " OUT SIG  090909009 "<<std::endl;
       }  
    }
    else {
      ntk.foreach_fanout(n, [this, &out_node](node const& nd){ out_node = nd;});
      auto children_out_node = ordered_children(out_node);
      std::cout << " out node : "<< children_out_node[0].n << std::endl;
      for(int i = 0; i < 1; i++){
      if(children_out_node[i].n == n && ntk.is_complemented(children_out_node[i].s)){
        //out_sig = children_out_node[i].s;
        out_sig_complemented = true;
        std::cout << " out node :INNN "<< children_out_node[0].n << std::endl;
      }
     }  
    }
    
    if(ntk.is_complemented(children[0].s) && ntk.is_complemented(children[1].s) && out_sig_complemented){
      return true;
    }
    return false;
  }
  // auto child_sig0 = ntk.create_and(children[1].s, grand_children0[1].s);
  //   auto grand_child_sig1 = ntk.create_and(children[1].s, great_grand_children0_0[1].s);
  //   ntk.substitute_node(great_grand_children0_0[1].n, grand_child_sig1);
  //   ntk.substitute_node(grand_children0[1].n, child_sig0);
  //   ntk.substitute_node(n, children[0].s);
  bool optimise_three_layer_distri(node &n, Node_sig &move_child, Node_sig &move_gc, Node_sig &move_ggc, Node_sig &main_child)
  {
    auto child_sig0 = ntk.create_and(move_child.s, move_gc.s);
    auto grand_child_sig1 = ntk.create_and(move_child.s, move_ggc.s);
    ntk.substitute_node(move_ggc.n, grand_child_sig1);
    ntk.substitute_node(move_gc.n, child_sig0);
    ntk.substitute_node(n, main_child.s);
    return true;
  }
  
  /* Try the distributivity rule on node n. Return true if the network is updated. */
  bool try_distributivity (node n)
 {
    /* TODO */
    std::cout << " NODE DISTRI: "<< n << std::endl; 
    
    auto children = ordered_children(n);
    bool both_pi = (ntk.is_pi(children[0].n)) & (ntk.is_pi(children[1].n));// nothing to do if both children are pi
      if(both_pi){
       std::cout << "##### NO UPDATE ******** BOTH PI " << std::endl;
       return false;
      }
    
    // both are nodes
    bool common_gc_exist;
    Node_sig common_gc;
    std::cout << " ^-^-^-^-^-^-^-^-^-^-^ CHECKING AIG OR::  "<< is_aig_or(n) << " "<< ntk.is_complemented(children[0].s)<< " "<<ntk.is_complemented(children[1].s) <<std::endl;
    //is_aig_or(n) && 
    bool both_child_complemented = (ntk.is_complemented(children[0].s)) & (ntk.is_complemented(children[1].s));
    
    if((is_aig_or(n) && ntk.is_on_critical_path(n)) || ( both_child_complemented && ntk.is_on_critical_path(n))){
      bool one_pi = (ntk.is_pi(children[0].n)) | (ntk.is_pi(children[1].n));// nothing to do if both children are pi
      if(one_pi){
       std::cout << "##### NO UPDATE ******** ALEAST ONE PI " << std::endl;
       return false;
      }
      std::cout << " HERE 1"<< std::endl;
      //if(ntk.is_and(children[0].n) & ntk.is_and(children[1].n)){
        auto grand_children0 = ordered_children(children[0].n);
        auto grand_children1 = ordered_children(children[1].n);
        bool all_gc_signal_same = (ntk.is_complemented(grand_children0[0].s) ^ ntk.is_complemented(grand_children0[1].s))
                                    ^ (ntk.is_complemented(grand_children1[0].s) ^ ntk.is_complemented(grand_children1[1].s));
        if(!all_gc_signal_same){
        // check if a common grand-child exists
        std::tie(common_gc, common_gc_exist) = check_common_grand_child(grand_children0, grand_children1);
        if(common_gc_exist == true && ntk.is_on_critical_path(common_gc.n)){
          
          auto non_common_gc = get_other_two_nodes(grand_children0, grand_children1, common_gc.n);
          
          auto hanging = ntk.create_and(ntk.create_not(non_common_gc[0].s), !non_common_gc[1].s);
          auto new_out = ntk.create_and(ntk.create_not(hanging), common_gc.s);
          std::cout<< "new_out signal : "<< ntk.is_complemented(new_out)<<std::endl;

          ntk.substitute_node(n, !new_out);
          std::cout << "~~~~~~~~~ UPDATE ~~~ common gc " << std::endl;
          return true;
        }
        else{
          std::cout << "##### NO UPDATE ******** NO COMMON GC " << std::endl;
          return false;
          

        }
      }
    }
    // three layer distributivity
    bool atlest_one_child_complemented = (ntk.is_complemented(children[0].s)) | (ntk.is_complemented(children[1].s));
    if( atlest_one_child_complemented && ntk.is_on_critical_path(n)){
      auto grand_children0 = ordered_children(children[0].n);
      auto grand_children1 = ordered_children(children[1].n);
      std::cout << "&&&& HERHER  1111" << std::endl;
      
      if(ntk.is_complemented(children[0].s) && ntk.is_on_critical_path(children[0].n)){
        bool both_gc0_signal_complemented = (ntk.is_complemented(grand_children0[0].s) & ntk.is_complemented(grand_children0[1].s));
        if(both_gc0_signal_complemented){
          auto great_grand_children0_0 = ordered_children(grand_children0[0].n);
          bool not_both_ggc0_signal_complemented = ~(ntk.is_complemented(great_grand_children0_0[0].s) & ntk.is_complemented(great_grand_children0_0[1].s));  
        
          auto great_grand_children0_1 = ordered_children(grand_children0[1].n);
          bool not_both_ggc1_signal_complemented = ~(ntk.is_complemented(great_grand_children0_1[0].s) & ntk.is_complemented(great_grand_children0_1[1].s));
          if(not_both_ggc0_signal_complemented && ntk.is_on_critical_path(grand_children0[0].n)) {

            Node_sig move_child = children[1];
            Node_sig move_gc = grand_children0[1];
            Node_sig move_ggc;
            Node_sig main_child = children[0];
            if(ntk.is_on_critical_path(great_grand_children0_0[0].n)){
              move_ggc = great_grand_children0_0[1];
            }
            else if (ntk.is_on_critical_path(great_grand_children0_0[1].n)){
              move_ggc = great_grand_children0_0[0];
            }
            
            return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child);
          }
          else if (not_both_ggc1_signal_complemented && ntk.is_on_critical_path(grand_children0[1].n)){
            Node_sig move_child = children[1];
            Node_sig move_gc = grand_children0[0];
            Node_sig move_ggc;
            Node_sig main_child = children[0];
            if(ntk.is_on_critical_path(great_grand_children0_1[0].n)){
              move_ggc = great_grand_children0_1[1];
            }
            else if (ntk.is_on_critical_path(great_grand_children0_1[1].n)){
              move_ggc = great_grand_children0_1[0];
            }
            return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child);
          }
          else {
            return false;
          }
        }
        else {
          return false;
        }
      }
      else if(ntk.is_complemented(children[1].s) && ntk.is_on_critical_path(children[1].n)){
        bool both_gc1_signal_complemented = (ntk.is_complemented(grand_children1[0].s) & ntk.is_complemented(grand_children1[1].s));
        if(both_gc1_signal_complemented){
          auto great_grand_children1_0 = ordered_children(grand_children1[0].n);
          bool not_both_ggc0_signal_complemented = ~(ntk.is_complemented(great_grand_children1_0[0].s) & ntk.is_complemented(great_grand_children1_0[1].s));  
        
          auto great_grand_children1_1 = ordered_children(grand_children1[1].n);
          bool not_both_ggc1_signal_complemented = ~(ntk.is_complemented(great_grand_children1_1[0].s) & ntk.is_complemented(great_grand_children1_1[1].s));
          if(not_both_ggc0_signal_complemented && ntk.is_on_critical_path(grand_children1[1].n)) {
            Node_sig move_child = children[0];
            Node_sig move_gc = grand_children1[0];
            Node_sig move_ggc;
            Node_sig main_child = children[1];
            if(ntk.is_on_critical_path(great_grand_children1_0[0].n)){
              move_ggc = great_grand_children1_0[1];
            }
            else if (ntk.is_on_critical_path(great_grand_children1_0[1].n)){
              move_ggc = great_grand_children1_0[0];
            }
            return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child);    
          }
          else if (not_both_ggc1_signal_complemented && ntk.is_on_critical_path(grand_children1[0].n)){
            Node_sig move_child = children[0];
            Node_sig move_gc = grand_children1[1];
            Node_sig move_ggc;
            Node_sig main_child = children[1];
            if(ntk.is_on_critical_path(great_grand_children1_1[0].n)){
              move_ggc = great_grand_children1_1[1];
            }
            else if (ntk.is_on_critical_path(great_grand_children1_1[1].n)){
              move_ggc = great_grand_children1_1[0];
            }
            return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child);  
          }
          else {
            return false;
          }
        
        }
        else {
          return false;
        }
      }
    }
  }

private:
  Ntk& ntk;
};

} // namespace detail

/* Entry point for users to call */
template<class Ntk>
void aig_algebraic_rewriting( Ntk& ntk )
{
  static_assert( std::is_same_v<typename Ntk::base_type, aig_network>, "Ntk is not an AIG" );

  depth_view dntk{ntk};
  fanout_view dntk2{dntk};
  detail::aig_algebraic_rewriting_impl p( dntk2 );
  p.run();
}

} /* namespace mockturtle */