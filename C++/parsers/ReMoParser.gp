header "pre_include_hpp"
{
    #include "bio/shared_parser_state.h"
    #include "bio/remo.h"
}

header "pre_include_cpp"
{
    #include "bio-pch.h"
    #include "bio/defs.h"
}

header "post_include_cpp"
{
    using namespace std;
    using namespace antlr;
    using namespace BIO_NS;
    
#ifdef _MSC_VER
        #pragma warning( disable : 4101 )
#endif
}

options
{
    language="Cpp";
    namespace = "bio";
}





class ReMoParser extends Parser;
options
{
    codeGenMakeSwitchThreshold = 999;
    importVocab = ReMoLexer;
}

{
public:
    ReMoSharedParserState sps;
    
}




remo_extraction returns [ ReMoExtraction::ptr_t ext ]
{
    ext.reset(new ReMoExtraction);
    ReMoSequenceGroup::ptr_t sg;
}
    :
    ( sg = sequence_group { ext->sequence_groups.push_back(sg); } )+
    ;





    
sequence_group returns [ ReMoSequenceGroup::ptr_t sg ]
{
    sg.reset(new ReMoSequenceGroup);
    ReMoSequence::ptr_t seq;
    ReMoBundle::ptr_t rb;
    double belief_threshold;
}
    :

    (
        belief_threshold = double_value
        NEW_LINE
    )?
    
    (
        seq = sequence
        {
            sg->sequences.push_back(seq);
        }
    )+
    SEQUENCELISTEND NEW_LINE
    
    (
        rb = remo_bundle
        {
            BOOST_ASSERT(rb->remos.end() != rb->remos.find(rb->centre_sequence));
            BOOST_ASSERT(! rb->remos[rb->centre_sequence].empty());
            ReMoSubSequence
                key(
                    rb->centre_sequence,
                    rb->remos[rb->centre_sequence].begin()->get()->range);
            sg->remo_bundles[key] = rb;
        } 
    )*
    REMOLISTEND NEW_LINE
    ;





sequence returns [ ReMoSequence::ptr_t seq ]
{
    seq.reset(new ReMoSequence);
    ReMoLocation loc = REMO_LOC_UNDEFINED;
    ReMoExon ex;
    ReMoSpecies sp;
    int length;
    ReMoSequenceId seq_id_details;
}
    :
    
    // The sequence id
    { sps.push_lexer("string"); }
    seq_id:STRING { seq->id = seq_id->getText(); } NEW_LINE
    { sps.pop_lexer(); }
    
    length = integer_value { seq->length = length; } NEW_LINE
    
    loc = location { seq->location = loc; } NEW_LINE
    
    sequence_position[seq] NEW_LINE

    sp = species { seq->species = sp; } NEW_LINE

    ( ex = exon { seq->exons.push_back(ex); } )* EXONLISTEND NEW_LINE
    
    ;


sequence_position [ ReMoSequence::ptr_t seq ]
{
    int position;
}
    :
    ( NONE ) =>
    (
        NONE
        { seq->has_position = false; }
    )
    |
    (
        position = integer_value
        { seq->has_position = true; seq->position = position; }
    )
    ;


remo_bundle returns [ ReMoBundle::ptr_t rb ]
{
    rb.reset(new ReMoBundle);
    ReMo::ptr_t re;
}
    :
    (
        { sps.push_lexer("string"); }
        seq_id:STRING NEW_LINE
        { sps.pop_lexer(); }
        {
            //if we haven't assigned the centre sequence yet, this (the first one) must be it
            if ("" == rb->centre_sequence)
            {
                rb->centre_sequence = seq_id->getText();
            }
        }

        (
            re = remo
            {
                //put the remo in the list for this sequence
                rb->remos[seq_id->getText()].push_back(re);
            }
        )+
    )*
    ENDREMO NEW_LINE
    ;


remo returns [ ReMo::ptr_t re ]
{
    re.reset(new ReMo);
    int conservation, repeat_ratio;
    double belief;
}
    :
    
    range [ re->range ]
    
    // the target sequence
    target [ re ]
    
    conservation = integer_value { re->conservation = conservation; } NEW_LINE
    repeat_ratio = integer_value { re->repeat_ratio = repeat_ratio; } NEW_LINE
    belief = double_value { re->belief = belief; } NEW_LINE
    
    { sps.push_lexer("string"); }
    masked_seq:STRING { re->masked_sequence = masked_seq->getText(); } NEW_LINE
    unmasked_seq:STRING { re->unmasked_sequence = unmasked_seq->getText(); } NEW_LINE
    { sps.pop_lexer(); }

    ;



    
target [ ReMo::ptr_t re ]
    :
    ( TARGETSELF ) =>
    (
        TARGETSELF NEW_LINE
        {
            re->target_range = re->range;
        }
    )
    |
    (
        range [ re->target_range ]
    )
    ;


subsequence [ ReMoSubSequence & ss ]
    :
    
    // The sequence id
    { sps.push_lexer("string"); }
    seq_id:STRING { ss.id = seq_id->getText(); } NEW_LINE
    { sps.pop_lexer(); }
    
    // where it is in the sequence
    range [ ss.range ]
    
    ;



    
location returns [ ReMoLocation loc ]
{
    loc = REMO_LOC_UNDEFINED;
}
    :
    UPSTREAM { loc = REMO_LOC_UPSTREAM; }
    | DOWNSTREAM { loc = REMO_LOC_DOWNSTREAM; }
    | GENEREGION { loc = REMO_LOC_GENEREGION; }
    | UNDEFINED { loc = REMO_LOC_UNDEFINED; }
    ;







species returns [ ReMoSpecies s ]
{
    s = REMO_SPECIES_UNKNOWN;
}
    :
    CHICK { s = REMO_SPECIES_CHICK; }
    | CHIMP { s = REMO_SPECIES_CHIMP; }
    | CIONA { s = REMO_SPECIES_CIONA; }
    | COW { s = REMO_SPECIES_COW; }
    | DOG { s = REMO_SPECIES_DOG; }
    | FLY { s = REMO_SPECIES_FLY; }
    | FUGU { s = REMO_SPECIES_FUGU; }
    | HUMAN { s = REMO_SPECIES_HUMAN; }
    | MOUSE { s = REMO_SPECIES_MOUSE; }
    | OPOSSUM { s = REMO_SPECIES_OPOSSUM; }
    | RAT { s = REMO_SPECIES_RAT; }
    | TETRAODON { s = REMO_SPECIES_TETRAODON; }
    | XENOPUS { s = REMO_SPECIES_XENOPUS; }
    | ZEBRAFISH { s = REMO_SPECIES_ZEBRAFISH; }
    ;
    
    

exon returns [ ReMoExon ex ]
    :
    ei:STRING {
        static const boost::regex re( "([A-Z]+)([0-9]+)" );
        boost::cmatch what; 
        const std::string text = ei->getText().c_str();
        const bool result = boost::regex_match( text.c_str(), what, re );
        ex.id.value = result ? what[1] : text;
        ex.id.number = result ? boost::lexical_cast< unsigned >( what[2] ) : 0;
    } NEW_LINE
    range [ ex.range ]
    ;
    


range [ ReMoRange & range ]
{
    int start, end;
}
    :
    start = integer_value NEW_LINE
    end = integer_value NEW_LINE
    {
        range.start = start;
        range.end = end;
    }
    ;


integer_value returns [ int i ]
{
    i = 0;
}
    :
    { sps.push_lexer("number"); }
    n:NUMBER
    {
        stringstream ss;
        ss << n->getText();
        ss >> i;
    }
    { sps.pop_lexer(); }
    ;
    
    
    
    
double_value returns [ double d ]
{
    d = 0;
}
    :
    { sps.push_lexer("number"); }
    n:NUMBER
    {
        stringstream ss;
        ss << n->getText();
        ss >> d;
    }
    { sps.pop_lexer(); }
    ;
    
