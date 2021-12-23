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
      std::cout << " ------ Level: " << ntk.level(n) << "         -------  " << "Node: "<< n << std::endl;
      ntk.foreach_fanin(n, [this]( auto const& f) { std::cout << "        " << ntk.get_node(f) << " ; is signal compl:"<< ntk.is_complemented(f)<<std::endl;
                                                                          });
    });
    ntk.foreach_po([&](signal const &f){ std::cout<< "po signal : "<< ntk.get_node(f) << "  " << ntk.is_complemented(f)<<std::endl;});
    ntk.foreach_pi([&](node const &f){ std::cout<< "pi node after update: "<< f <<  std::endl;});
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
      std::cout << " ------ Level: " << ntk.level(n) << "         -------  " << "Node: "<< n << std::endl;
      ntk.foreach_fanin(n, [this]( auto const& f) { std::cout << "        " << ntk.get_node(f) << " ; is signal compl :"<< ntk.is_complemented(f)<< std::endl;
                                                                          });
    });
    ntk.foreach_po([&](signal const &f){ std::cout<< "po signal after update: "<< ntk.get_node(f) << "  " <<ntk.is_complemented(f)<<std::endl;});
    ntk.foreach_pi([&](node const &f){ std::cout<< "pi node after update: "<< f << std::endl;});

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
    bool out_sig_comp(node n)
  {
    auto children = ordered_children(n);
    node out_node;
    signal out_sig;
    bool po_node = false;
    bool out_sig_complemented = false;
       ntk.foreach_po([this, &po_node, &n, &out_sig, &out_sig_complemented](signal const& ns){ out_sig = ns;
       //std::cout<< " ^^^^^^^ HERE FOR PO  "<< ntk.get_node(out_sig) <<"  ==  "<< n<<  std::endl;
       if(ntk.get_node(out_sig) == n ){
         po_node = true;
         if(ntk.is_complemented(out_sig) == true){
         out_sig_complemented = true;
         //std::cout<< " ^^^^^^^ HERE FOR PO TREU   "<< ntk.get_node(out_sig) <<"  ==  "<< n<<  std::endl;
         return true;
         }
       }
       });
         
      // cannot use this when node does not have any parent node, in other words,  when node is connected to primary output
      if(po_node == false){
        ntk.foreach_fanout(n, [this, &out_node](node const& nd){ out_node = nd;});
        auto children_out_node = ordered_children(out_node);
        for(int i = 0; i < children_out_node.size(); i++){
          if(children_out_node[i].n == n && ntk.is_complemented(children_out_node[i].s)){
            out_sig_complemented = true;
            return true;
          }
        }
      } 
    return false;
  }
  std::tuple<Node_sig, bool, bool> check_common_grand_child(std::array<Node_sig, 2> grand_children0, std::array<Node_sig, 2> grand_children1)
  {
    bool polarity_check = true;
    bool atleast_xor = false;
    for(int gr0 = 0; gr0 < 2; gr0++){
      for(int gr1 = 0; gr1 < 2; gr1++){
        if(grand_children0[gr0].n == grand_children1[gr1].n){
          polarity_check = ntk.is_complemented(grand_children0[gr0].s) ^ ntk.is_complemented(grand_children1[gr1].s);
          atleast_xor = true;
          if(!polarity_check){
            return std::make_tuple(grand_children0[gr0], true, atleast_xor);
          }
        }
      }
    }
    return std::make_tuple(grand_children0[0], false, atleast_xor);
  }
  std::tuple<Node_sig, bool, bool> check_common_great_grand_child(std::array<Node_sig, 2> grand_children0, std::array<Node_sig, 2> grand_children1)
  {
    bool polarity_check = true;
    bool atleast_xor = false;
    for(int gr0 = 0; gr0 < 2; gr0++){
      for(int gr1 = 0; gr1 < 2; gr1++){
        if(grand_children0[gr0].n == grand_children1[gr1].n){
          polarity_check = ntk.is_complemented(grand_children0[gr0].s) ^ ntk.is_complemented(grand_children1[gr1].s);
          atleast_xor = true;
          if(!polarity_check){
            return std::make_tuple(grand_children0[gr0], true, atleast_xor);
          }
        }
      }
    }
    return std::make_tuple(grand_children0[0], false, atleast_xor);
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
    std::cout << " ** Entering Assosciativity for Node : "<< n << std::endl; 
    Node_sig c_on_critical_path, c_not_on_critical_path;
    Node_sig gc_on_critical_path, gc_not_on_critical_path;
    bool com_gc, xor_op;
    Node_sig hi;
    auto children = ordered_children(n);
    bool only_one_child_comp = (ntk.is_complemented(children[0].s)) ^ (ntk.is_complemented(children[1].s));
    
    // for and associativity
    bool one_pi = (ntk.is_pi(children[0].n)) & (ntk.is_pi(children[1].n));//nothing to do if both children are pi or both nodes
    if(one_pi){
      //std::cout << "          NO UPDATE: BOTH PI " << std::endl;
      return false;
    } 
    if(ntk.is_on_critical_path(n)){
      if(!ntk.is_pi(children[0].n) && children[0].cp == true && children[1].cp == false){
        c_on_critical_path = children[0];
        c_not_on_critical_path = children[1];
      }
      else if(!ntk.is_pi(children[1].n) && children[1].cp == true && children[0].cp == false){
        c_on_critical_path = children[1];
        c_not_on_critical_path = children[0];
      }
      else if(!ntk.is_pi(children[0].n) && !ntk.is_pi(children[1].n) && children[0].cp == true && children[1].cp == true){
        return false;
      }
      
      else {
        //std::cout << "          NO UPDATE: No child on CP " << std::endl;
        return false;
      }
      
      if(!ntk.is_complemented(c_on_critical_path.s)){
        auto grand_children = ordered_children(c_on_critical_path.n);
        auto grand_children_not_cp = ordered_children(c_not_on_critical_path.n);
        
        if(grand_children[0].cp == true && grand_children[1].cp == false){
          gc_on_critical_path = grand_children[0];
          gc_not_on_critical_path = grand_children[1];
        }
        else if(grand_children[1].cp == true && grand_children[0].cp == false){
          gc_on_critical_path = grand_children[1];
          gc_not_on_critical_path = grand_children[0];
        }
        else if(grand_children[0].cp == true && grand_children[1].cp == true){
          return false;
        }
      
        else {
          //std::cout << "          NO UPDATE: No grand-child on CP " << std::endl;
          return false;
        }

        //std::tie(std::ignore, com_gc, std::ignore) = check_common_grand_child(grand_children, grand_children_not_cp);
        //std::tie(std::ignore, std::ignore, xor_op) = check_common_great_grand_child(ordered_children(grand_children[0].n), ordered_children(grand_children[1].n));
        if((ntk.level(gc_on_critical_path.n) >= ntk.level(c_not_on_critical_path.n) + 1) ){//&& !xor_op
          auto hanging = ntk.create_and(gc_not_on_critical_path.s, c_not_on_critical_path.s);
          auto new_out = ntk.create_and(gc_on_critical_path.s, hanging);
          ntk.substitute_node(n, new_out);
          std::cout << "- - - - - - UPDATE: OR Associativity " << "  COMMON GC "<< com_gc << std::endl;
          return true;
        }
        else{
          //std::cout << "          NO UPDATE: Level difference less than or equal to 1 " << std::endl;
          return false;
        }
      }
      else 
        return false;

    }
    
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


  bool optimise_three_layer_distri(node &n, Node_sig &move_child, Node_sig &move_gc, Node_sig &move_ggc, Node_sig &main_child, Node_sig &main_gchild, Node_sig &main_ggc)
  {
    std::cout << " ======== OPTIMISE "<< move_child.n << " " << move_gc.n << std::endl;
    std::cout << " ======== OPTIMISE "<< move_child.n << " " << move_ggc.n << std::endl;
    std::cout << " ======== OPTIMISE IS GC COMP : "<< ntk.is_complemented(move_gc.s) << std::endl;
    
    auto grand_child_sig1 = ntk.create_and(move_child.s, move_ggc.s);
    auto child_sig1 = ntk.create_and(grand_child_sig1, main_ggc.s);

    auto child_sig0 = ntk.create_and(move_child.s, !move_gc.s);
    
    auto new_out = ntk.create_and(!child_sig0, !child_sig1);
    std::cout << " ======== OPTIMISE  SUB " << move_ggc.n << std::endl;
    std::cout << " ======== OPTIMISE  SUB " << move_gc.n << std::endl;
    ntk.substitute_node(n, !new_out);
    return true;
  }
  
  /* Try the distributivity rule on node n. Return true if the network is updated. */
  bool try_distributivity (node n)
 {
    /* TODO */
    std::cout << " ++ Entering Disttributivity for Node : "<< n << std::endl; 
    
    auto children = ordered_children(n);
    bool both_pi = (ntk.is_pi(children[0].n)) & (ntk.is_pi(children[1].n));// nothing to do if both children are pi
    if(both_pi){
       //std::cout << "          NO UPDATE: BOTH PI " << std::endl;
       return false;
    }
    
    // both are nodes
    bool common_gc_exist;
    Node_sig common_gc;

    bool both_child_complemented = (ntk.is_complemented(children[0].s)) & (ntk.is_complemented(children[1].s));
    
    bool atlest_one_child_complemented = (ntk.is_complemented(children[0].s)) | (ntk.is_complemented(children[1].s));
    Node_sig move_ggc;
    Node_sig ggc_on_cp;
    Node_sig ggc_not_cp;

    if((out_sig_comp(n) && both_child_complemented && ntk.is_on_critical_path(n)) || (!out_sig_comp(n) && both_child_complemented && ntk.is_on_critical_path(n) )){
      bool one_pi = (ntk.is_pi(children[0].n)) | (ntk.is_pi(children[1].n));// nothing to do if both children are pi
      if(one_pi){
       //std::cout << "          NO UPDATE: AND ON TOP/ OR ON TOP, There is ATLEAST ONE PI " << std::endl;
       return false;
      }
      //std::cout << " HERE 1"<< std::endl;
      //if(ntk.is_and(children[0].n) & ntk.is_and(children[1].n)){
        auto grand_children0 = ordered_children(children[0].n);
        auto grand_children1 = ordered_children(children[1].n);
        bool all_gc_signal_same = (ntk.is_complemented(grand_children0[0].s) ^ ntk.is_complemented(grand_children0[1].s))
                                    ^ (ntk.is_complemented(grand_children1[0].s) ^ ntk.is_complemented(grand_children1[1].s));
        //if(!all_gc_signal_same){
        // check if a common grand-child exists
        std::tie(common_gc, common_gc_exist, std::ignore) = check_common_grand_child(grand_children0, grand_children1);
        if(common_gc_exist == true && ntk.is_on_critical_path(common_gc.n)){
          
          auto non_common_gc = get_other_two_nodes(grand_children0, grand_children1, common_gc.n);
          
          auto hanging = ntk.create_and(!non_common_gc[0].s, !non_common_gc[1].s);
          auto new_out = ntk.create_and(!hanging, common_gc.s);
          //std::cout<< "new_out signal : "<< ntk.is_complemented(new_out)<<std::endl;

          ntk.substitute_node(n, !new_out);
          std::cout << "- - - - - UPDATE: COMMON grand-child found " << std::endl;
          return true;
        }
        else{
          //std::cout << "          NO UPDATE: No common grand-child " << std::endl;
          return false;
          

        }
      }
    //}
    // three layer distributivity
// && atlest_one_child_complemented 
    else if( !out_sig_comp(n)&& ntk.is_on_critical_path(n)){
      auto grand_children0 = ordered_children(children[0].n);
      auto grand_children1 = ordered_children(children[1].n);
      Node_sig dc;
      bool common_ggc;
      if(ntk.is_on_critical_path(children[0].n) && ntk.is_on_critical_path(children[1].n)){
        return false;
      }

      if(ntk.is_complemented(children[0].s) && ntk.is_on_critical_path(children[0].n)  && !ntk.is_pi(children[1].n)){
        Node_sig child_on_cp = children[0];
        Node_sig child_not_cp = children[1];
        bool both_gc0_signal_complemented = (ntk.is_complemented(grand_children0[0].s) & ntk.is_complemented(grand_children0[1].s));
        if(both_gc0_signal_complemented){
          auto great_grand_children0_0 = ordered_children(grand_children0[0].n);
          //bool both_ggc0_signal_complemented = ntk.is_complemented(great_grand_children0_0[0].s) & ntk.is_complemented(great_grand_children0_0[1].s);  
        
          auto great_grand_children0_1 = ordered_children(grand_children0[1].n);
          if(ntk.is_on_critical_path(grand_children0[0].n) && ntk.is_on_critical_path(grand_children0[1].n)){
            return false;
          }
          //std::tie(dc, common_ggc) = check_common_grand_child(great_grand_children0_0, great_grand_children0_1);
          //bool both_ggc1_signal_complemented = ntk.is_complemented(great_grand_children0_1[0].s) & ntk.is_complemented(great_grand_children0_1[1].s);
          if(ntk.is_on_critical_path(grand_children0[0].n) && !ntk.is_pi(grand_children0[0].n)) {
            Node_sig gc_on_cp = grand_children0[0];
            Node_sig gc_not_cp = grand_children0[1];

            Node_sig move_child = child_not_cp;
            Node_sig move_gc = gc_not_cp;
            
            Node_sig main_child = child_on_cp;
            //  && !ntk.is_pi(great_grand_children0_0[0].n)
            if(ntk.is_on_critical_path(great_grand_children0_0[0].n) && !ntk.is_on_critical_path(great_grand_children0_0[1].n)){
              ggc_on_cp = great_grand_children0_0[0];
              ggc_not_cp = great_grand_children0_0[1];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 1 " << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);
            }// && !ntk.is_pi(great_grand_children0_0[1].n)
            else if (ntk.is_on_critical_path(great_grand_children0_0[1].n) && !ntk.is_on_critical_path(great_grand_children0_0[0].n)){
              ggc_on_cp = great_grand_children0_0[1];
              ggc_not_cp = great_grand_children0_0[0];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 2" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);
            }
            else {
              //std::cout << "          NO UPDATE: Three layer distributivity - 3" << std::endl;
              return false;
            }
            
            
          }
          else if (ntk.is_on_critical_path(grand_children0[1].n) && !ntk.is_pi(grand_children0[1].n)){
            Node_sig gc_on_cp = grand_children0[1];
            Node_sig gc_not_cp = grand_children0[0];

            Node_sig move_child = child_not_cp;
            Node_sig move_gc = gc_not_cp;
            
            Node_sig main_child = child_on_cp;// && !ntk.is_pi(great_grand_children0_1[0].n)
            if(ntk.is_on_critical_path(great_grand_children0_1[0].n) && !ntk.is_on_critical_path(great_grand_children0_1[1].n)){
              ggc_on_cp = great_grand_children0_1[0];
              ggc_not_cp = great_grand_children0_1[1];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 4" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);
            }// && !ntk.is_pi(great_grand_children0_1[1].n)
            else if (ntk.is_on_critical_path(great_grand_children0_1[1].n) && !ntk.is_on_critical_path(great_grand_children0_1[0].n)){
              ggc_on_cp = great_grand_children0_1[1];
              ggc_not_cp = great_grand_children0_1[0];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 5" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);
            }
            else {
              //std::cout << "          NO UPDATE: Three layer distributivity - 6" << std::endl;
              return false;
            }
            
            
          }
          else {
            //std::cout << "          NO UPDATE: Three layer distributivity - 7" << std::endl;
            return false;
          }
        }
        else {
          //std::cout << "          NO UPDATE: Three layer distributivity - 8" << std::endl;
          return false;
        }
      }
      else if(ntk.is_complemented(children[1].s) && ntk.is_on_critical_path(children[1].n)  && !ntk.is_pi(children[1].n)){
        bool both_gc1_signal_complemented = (ntk.is_complemented(grand_children1[0].s) & ntk.is_complemented(grand_children1[1].s));
        Node_sig child_on_cp = children[1];
        Node_sig child_not_cp = children[0];
        
        if(both_gc1_signal_complemented){
          auto great_grand_children1_0 = ordered_children(grand_children1[0].n);
          bool both_ggc0_signal_complemented = ntk.is_complemented(great_grand_children1_0[0].s) & ntk.is_complemented(great_grand_children1_0[1].s);  
          auto great_grand_children1_1 = ordered_children(grand_children1[1].n);
          //std::tie(dc, common_ggc) = check_common_grand_child(great_grand_children1_0, great_grand_children1_1);
          bool both_ggc1_signal_complemented = ntk.is_complemented(great_grand_children1_1[0].s) & ntk.is_complemented(great_grand_children1_1[1].s);
          if(ntk.is_on_critical_path(grand_children1[0].n) && ntk.is_on_critical_path(grand_children1[1].n)){
            return false;
          }
          if(ntk.is_on_critical_path(grand_children1[0].n) && !ntk.is_pi(grand_children1[0].n)) {
            Node_sig gc_on_cp = grand_children1[0];
            Node_sig gc_not_cp = grand_children1[1];

            Node_sig move_child = child_not_cp;
            Node_sig move_gc = gc_not_cp;
            Node_sig main_child = child_on_cp;// && !ntk.is_pi(great_grand_children1_0[0].n)
            if(ntk.is_on_critical_path(great_grand_children1_0[0].n) && !ntk.is_on_critical_path(great_grand_children1_0[1].n)){
              ggc_on_cp = great_grand_children1_0[0];
              ggc_not_cp = great_grand_children1_0[1];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 9" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);    
            }//  && !ntk.is_pi(great_grand_children1_0[1].n)
            else if (ntk.is_on_critical_path(great_grand_children1_0[1].n) && !ntk.is_on_critical_path(great_grand_children1_0[0].n)){
              ggc_on_cp = great_grand_children1_0[1];
              ggc_not_cp = great_grand_children1_0[0];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 10" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);    
            }
            else {
              //std::cout << "          NO UPDATE: Three layer distributivity - 11" << std::endl;
              return false;
            }
            
          }
          else if (ntk.is_on_critical_path(grand_children1[1].n)  && !ntk.is_pi(grand_children1[1].n)){
            //auto great_grand_children1_1 = ordered_children(grand_children1[1].n);
            Node_sig gc_on_cp = grand_children1[1];
            Node_sig gc_not_cp = grand_children1[0];

            Node_sig move_child = child_not_cp;
            Node_sig move_gc = gc_not_cp;
            
            Node_sig main_child = child_on_cp;
            // && !ntk.is_pi(great_grand_children1_1[0].n)
            if(ntk.is_on_critical_path(great_grand_children1_1[0].n) && !ntk.is_on_critical_path(great_grand_children1_1[1].n)){
              ggc_on_cp = great_grand_children1_1[0];
              ggc_not_cp = great_grand_children1_1[1];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 12" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);  
            }// && !ntk.is_pi(great_grand_children1_1[1].n)
            else if (ntk.is_on_critical_path(great_grand_children1_1[1].n) && !ntk.is_on_critical_path(great_grand_children1_1[0].n) ){
              ggc_on_cp = great_grand_children1_1[1];
              ggc_not_cp = great_grand_children1_1[0];
              move_ggc = ggc_not_cp;
              std::cout << "- - - - - UPDATE: Three layer distributivity - 13" << std::endl;
              return optimise_three_layer_distri(n, move_child, move_gc, move_ggc, main_child, gc_on_cp, ggc_on_cp);  
            }
            else {
              //std::cout << "          NO UPDATE: Three layer distributivity - 14" << std::endl;
              return false;
            }
          }
          else {
            //std::cout << "          NO UPDATE: Three layer distributivity - 15" << std::endl;
            return false;
          }
        
        }
        else {
          //std::cout << "          NO UPDATE: Three layer distributivity - 16" << std::endl;
          return false;
        }
      }
    }
    else {
      //std::cout << "          NO UPDATE: For distributivity " << std::endl;
      return false;
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