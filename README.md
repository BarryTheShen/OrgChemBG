# Organic Chemistry Board Game

[中文](https://github.com/BarryTheShen/OrgChemBG/blob/main/README_CHI.md) | English

**Special Thanks**  
Special thanks to [nandeck](http://www.nandeck.com/), Wikipedia, and ChatGPT for their invaluable contributions in composing the cards, providing some of the reactant images, and assisting with coding and game development.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. You are free to share and adapt the material as long as you give appropriate credit, provide a link to the license, and indicate if changes were made. You may not use the material for commercial purposes.

## Introduction

This is a board game for organic chemistry **enthusiasts**. This entire board game is inspired by the idea of total synthesis and how it is similar to the idea of traditional puzzles where the problem that needs to be solved is the process to make the final product. However, due to the complex nature of organic chemistry total synthesis, we as the creator team decided to not consider the following to make this board game playable and simple enough so that it can be done in a short period of time:

- Stereochemistry
- Solvents
- Temperature of the reactions (except for heating)
- Aromaticity (added in the future)
- Regioselectivity (If there are controversies, the player that does the synthesis decides)
- Desired yield

To make this playable, some places is made better for playability and many scientific rigidity is reduced.

## Objective of the Game

The objective of the game is the first to synthesize the final product that is given to you.

There are a few game modes that can be implemented

- **Non-Competitive:** This game mode is where everyone gets 3 different product and the first person to synthesize one of their products will win. Everyone else have the rest of the round to complete their product if they can.
- **Competitive:** This game mode consists of 3 final products that is shared by everyone. The first person to synthesize one of the product will win the game.
- **Point Based (Not implemented):** This game mode is currently not implemented and in the future players will have the chance to synthesize the most amount of products and each product will be assigned points to them and the person with the most amount of points will win the game.

## Components

There are a total of 4 classes of cards: Reaction cards, PDT, Premium, and Final Products

### Reaction cards

The reaction cards are responsible for synthesizing anything in the board game. They are cards that can be used as the “left side” of the reaction. You can use this to modify a current intermediate you have or you can initiate a new reaction. There are many subclasses of reaction cards:

- Reactant (dark yellow): These cards are used as your reactant of the reaction. These will be used to be used with the other reagents to synthesize new materials, and the newly synthesize materials can be modified with other reactants or reagents.
- Oxidant (dark red): This set of cards will be primarily be responsible for oxidizing your reactants. However this doesn’t include any acids that might be responsible for any oxidative cleavage or anything similar.
- Reductant (pink): This set of cards will be primarily be responsible for reducing your reactants. However, reductants that are metals are not included in the reductant category but rather in the Metals & Ylides section.
- Acid (light green): This set off cards are all acids in the deck including all the acids that might be responsible for oxidative cleavage or halogenation. It doesn’t include Lewis acids, but only Bronsted-Lowry acids.
- Base (cyan): This set of cards will be all the base in the deck.
- Halogen (purple): This set of cards will be responsible for modifying the reactant with halogens. However, this doesn’t include any acids (HBr) that might also be responsible for halogenation.
- Metals & Ylides (silver): This set of cards are a combination of metals and ylides. Some metals that are reductants or catalyzes other reactions are put in this category.
- Protecting & Leaving groups (Represented as PROT & LEAV GROUPS in game, dark green): This includes molecules that are either groups that can protect a function group or a modification to better the leaving group.
- Reaction condition & H2O (Represented as RXN COND & H2O, orange): This includes reaction condition cards (heating up, radical initiators) and also water.
- Other (dark grey): This includes a wide range of products that might not fit into any of the previous categories or there are too less to start a new pile. There are dehydration reagents, some special reaction compounds, and some random reagents.
- Note: All of these cards can be used interchangeably. A card that is in the Protecting & Leaving group section can be used as a reactant. They are not restricted to be a protecting group and can be used as a reactant. The categorization in this case is only for selecting the cards when it is a player’s turn.

### PDT

PDT is a class of cards that are small products that gives you access to premium class of cards when you synthesize them. They act like a currency. Whenever you synthesize a PDT, you get that card and you take it from the section. PDTs are designed to be synthesized in a 1 step 2-3 reactant/reagent reaction. This means that with the right material you can get the PDT a single step. There are a total of 48 PDT cards that can be synthesized.

### Premium

These are the premium cards or actions that can be bought with PDTs and used to assist your synthesis process:

- Wild Card (cost 1 PDT): Players can draw any reactant that is less than 2 carbon or any reagent that adds less than 2 carbon. Players can also draw any reactant or reagent that is already in the deck.
- Duplication Card (cost 1 PDT): Players can duplicate any thing they have and make it 2 cards. This includes their current molecule that they have or any other reagent that they want to duplicate.
- Reserve (cost 2 PDT)(not implemented): Players can reserve a final product that they want to synthesize.

### Final Products

These are the ultimate end goal of the game where players need to synthesize them to either win the game or get points.

### Extras

Other than these cards that is provided in the deck. some pens and papers is also recommended as a player need to draw out. This can be substituted with chemical structures to represent the chemicals synthesized.

## Setup

1. Separate the cards according to their color. For all the reaction cards, place them on top of the table side by side. The PDT, Premium, and Final Product cards should be placed on the right side of the desk. All cards should be facing downwards and separated into their designated sections.
2. For the PDT cards, pick 5 up and place them on the table facing upwards.
3. For Premium cards, flip them over and separate them into 2 piles: a Wild card pile and a Duplication card pile.
4. (Only for **Competitive game mode**): pick up 4 final product cards, and place them on the table facing upwards.
5. Each person should have a pen and a piece of paper in front of them that they use as their “synthesis paper”. This is you will use it to draw out your products and each time you do an reaction you are just redrawing this to change it. However, this can be substituted for a chemistry model, but scratch paper is suggested as total synthesis also often comes with lots of retrosynthesis.

This is an example of what it should look like for the non-competitive game mode:

![ Alt Text](https://github.com/BarryTheShen/OrgChemBG/blob/main/Resources/noncompetitive%20image.png)

This is an example of what it should look like for the competitive game mode:

![ Alt Text](https://github.com/BarryTheShen/OrgChemBG/blob/main/Resources/competitive%20image.png)

## Gameplay

The game is played over multiple rounds each person needs to go through the following phases to move on to the next person:

1. Action
2. Discard
3. Clean

### Action

This is where the player can do something. The player can choose one of the following actions:

- Draw
- React

**Draw**

- This is where the player can gather their reagents for their future reactions to make their reactions come true.
- The player can choose the following:
  1. Take 3 cards from 3 different reaction cards classes
  2. Take 2 cards from the same pile of reaction cards.
  3. *Steal 1 card from a player (draw a random card from them) (PDTs and Premium cards cannot be stolen)

*Note: For every round, this can be only done once for each player. If Peter stole 1 card from Andrew, then Josh can’t steal 1 card from Andrew within 1 round after Peter stole Andrew’s card. Peter can only steal a card from Andrew 2 rounds after.

If there are no cards left from a pile, the player need to choose another pile to take from.

- After drawing the card, the player moves on to the next phase: discard.

**React**

- This is where the fun of total synthesis really comes in. With the player’s current knowledge of chemistry, they have to try to make a valid reaction with their reactants. Each card represents 1 equivalent of the reactant/reagent that is on there. For example, if a reaction requires, 2 equivalent of a reactant/reagent, the player can get 2 of that card or use a wild card or duplication card instead.
- The player will place the card down that they want to react in the center of the table and explain their reaction. If there are problems with regioselectivity, the person who played the card will have the final say in the product, unless it is very obvious that reaction will not happen because of other reactions that will definitely take priority.
- After the reaction is completed, the player will draw their product on their paper. The player can then choose to redeem the product or not (declare a win or claim a PDT). If the player produces 2 PDTs at the same time they can claim them both.
- Although PDTs are usually produced at the end of React, PDTs can be redeemed any time within the player’s turn. Premium cards apply the same rule.
- If the reaction involves some sort of cleavage (ozonolysis) the player can choose to keep both of the products.
- After the reaction and redeeming is completed, proceed to the next phase: discard.

### Discard

- If a player has more than 7 Reaction cards (excluding Premium cards and PDTs), they need to discard the cards until they have less than or equal to 7 cards.

### Clean

- This is where the board is reset to ensure the next person can proceed without a problem.
- All the cards that is played and/or discarded is placed back to their designated piles and reshuffled so that the reagents are replenished.
- For PDT cards that’s claimed, replenish their spots with new PDT cards drawn from the pile, so the board looks like the starting position.
- After cleaning and resetting the board, it is the next person’s turn.

## Game’s End

- The game ends with one person successfully synthesizing the product. The rest of the round can be finished. (If Peter started first, then Andrew, and then Josh, and Andrew finishes the product. Josh will have 1 more turn and then the game ends. Peter doesn’t get to go again since he started first.)
- If there is a tie between 2 people. The person with the most synthesized products will win (if someone synthesize 2 products at the same time). If both players synthesize the same amount of products, then the person with the most PDTs and Premium cards win. If both players have the same amount for PDTs and Premium cards, the person with most amount of Reaction cards wins.

## Player Aids

### Category of cards and amount

| Reactant             | Oxidant      | Reductant              | Acid            | Base            |
| -------------------- | ------------ | ---------------------- | --------------- | --------------- |
| Ethanol (3)          | Ag₂O (1)     | B₂H₆ (1)               | Acetic acid (3) | DBU (1)         |
| Methanol (3)         | CrO₃ (2)     | CeCl₃ (1)              | H₂SO₄ (3)       | LDA (2)         |
| Glycerol (3)         | DMP (1)      | DIBAL (4)              | H₃PO₄ (1)       | NaNH₂ (2)       |
| But-1,3-diene (3)    | H₂CrO₄ (1)   | H₂ (6)                 | HBr (2)         | NaOH (6)        |
| Maleic Anhydride (3) | H₂O₂ (3)     | LiAlH(OtBu)₃ (1)       | HCl (5)         | n-BuLi (2)      |
| Urea (3)             | Hg(OAc)₂ (2) | LiAlH₄ (9)             | HI (2)          | Pyridine (2)    |
| Chloromethane (2)    | HgSO₄ (1)    | Lindler's catalyst (1) | HIO₄ (1)        | Pyrrolidine (1) |
| Bromomethane (2)     | KMnO₄ (7)    | NaBH₄ (2)              | HNO₃ (2)        | NaH (1)         |
| Iodomethane (2)      | mCPBA (2)    | Ni₂B (1)               | Oxalic acid (3) | K₂CO₃ (1)       |
| Acetone (3)          | MnO₂ (1)     |                        | TFA (1)         |                 |
| Propene (3)          | NaClO₂ (1)   |                        |                 |                 |
| Nitromethane (2)     | NaIO₄ (1)    |                        |                 |                 |
| Acetyl chloride (3)  | O₃ (3)       |                        |                 |                 |
| Propyne (3)          | OsO₄ (1)     |                        |                 |                 |
|                      | Pb(OAc)₄ (1) |                        |                 |                 |
|                      | PCC (1)      |                        |                 |                 |
|                      | SeO₂ (1)     |                        |                 |                 |

| Halogen   | Metals & Ylides          | Protecting & Leaving groups | Reaction condition & H₂O | Other     |
| --------- | ------------------------ | --------------------------- | ------------------------ | --------- |
| Br₂ (4)   | Li (3)                   | BocCl (1)                   | ∆ (4)                    | BF₃ (1)   |
| Cl₂ (2)   | Mg (3)                   | BsCl (1)                    | H₂O (7)                  | CH₂N₂ (1) |
| I₂ (2)    | Na (2)                   | Ethylene glycol (3)         | hv (3)                   | CO (2)    |
| KI (1)    | Ni (1)                   | MsCl (1)                    | AIBN (2)                 | CO₂ (2)   |
| NBS (1)   | Pd/C (1)                 | NaN₃ (1)                    |                          | DCC (3)   |
| NCS (1)   | PPh₃ (3)                 | TMSCl (1)                   |                          | EDT (1)   |
| NIS (1)   | Pt (1)                   | TsCl (1)                    |                          | KCN (3)   |
| PBr₃ (1)  | Ra-Ni (2)                |                             | N₂H₄ (1)                 |           |
| PCl₃ (1)  | Sulfur ylide (1)         |                             | NH₂OH (1)                |           |
| PCl₅ (1)  | Zn (Or any Zn alloy) (3) |                             | P₂O₅ (1)                 |           |
| SOBr₂ (1) |                          |                             | POCl₃ (1)                |           |
| SOCl₂ (1) |                          |                             |                          |           |
